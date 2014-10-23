suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

# bamtoolsr.r - some commands for processing BAM files

option_list <- list(
  make_option(c("-f", "--file"), 
              type="character",
              default=NA,
              help="Path of BAM file to process"),
  make_option(c("-e", "--extension"),
              type="character",
              default="native",
              help="Extension length ('native' for no extension, 'auto' to estimate fragment size)"),
  make_option(c("-n", "--name"),
              type="character",
              help="R variable name for resulting coverage and ranges objects"),
  make_option(c("-b", "--bigwig"),
              type="character",
              default=NA,
              help="File name for resulting BigWig"),
  make_option(c("-p", "--paired"),
              action="store_true",
              default=FALSE,
              help="Paired-end data"),
  make_option(c("-s", "--skipfilter"),
              action="store_true",
              default=FALSE,
              help="Skip artifact filter")
  )

readBAM <- function(filepath) {
  bam.gr <- as(readGAlignments(filepath), "GRanges")
  bam.gr
}

readBAM_pe <- function(filepath) {
  bam.gr <- granges(readGAlignmentPairs(filepath))
  strand(bam.gr) <- "*"
  bam.gr
}

write_bigwig <- function(cov, filename) {
  message(id, "Writing bigWig...")
  export(cov, filename)
}

artifact_filter <- function(bam.gr, dup_count_limit=10, ext_length) {
  reads <- data.table(as.data.frame(bam.gr))
  setnames(reads, c("chr", "start", "end", "width", "strand"))
  setkey(reads, chr, start, end, strand)
    
  reads$rownum <- 1:nrow(reads)

  filtered.reads <- NULL
  message(id, "artifact_filter: Starting with ", pn(nrow(reads)), " reads")
  for(chromosome in sort(unique(reads$chr))) {
    nothing <- gc()
    
    reads.chr     <- reads[chr == chromosome]
    reads.pos.chr <- reads.chr[strand == "+"]
    reads.neg.chr <- reads.chr[strand == "-"]
      
    ranges.pos.chr <- IRanges(start=reads.pos.chr$start, end=reads.pos.chr$start)
    ranges.neg.chr <- IRanges(start=reads.neg.chr$end,   end=reads.neg.chr$end)
    
    reads.dc <- reads.chr[, list(dupcount=length(rownum)), by=list(chr, start, end, strand)]

    dup.targets         <- reads.dc[dupcount >  dup_count_limit]
    dup.below_threshold <- reads.dc[dupcount <= dup_count_limit]
    
    if(nrow(dup.targets) == 0) {
      message(id, "(", chromosome, ") No duplicate reads exceed threshold of ", dup_count_limit)
      filtered.reads <- rbind(filtered.reads, reads.chr[, list(chr, start, end, strand)])
      next
    }
 
    # return all reads below duplicated threshold
    dup.below_threshold.clean  <- dup.below_threshold[, list(chr, start, end, strand)]
    dup.below_threshold.counts <- dup.below_threshold$dupcount 

    filtered.reads <- rbind(filtered.reads, dup.below_threshold.clean[rep(1:nrow(dup.below_threshold.clean), times=dup.below_threshold.counts)])
    rm(dup.below_threshold.clean)
    rm(dup.below_threshold)
    nothing <- gc()
    
    dup.targets$paired.start <- dup.targets[, ifelse(strand == "+", start + ext_length - 25, end - ext_length - 25)]
    dup.targets$paired.end   <- dup.targets[, ifelse(strand == "+", start + ext_length + 25, end - ext_length + 25)]
    
    pos_dups <- which(dup.targets$strand == "+")
    neg_dups <- which(dup.targets$strand == "-")
    
    dup.targets$paired.region.counts <- 0
    
    if(length(pos_dups) > 0) {
      dup.targets[pos_dups, ]$paired.region.counts <- counts_for_ranges(dup.targets[pos_dups, ], "paired.start", "paired.end", ranges.neg.chr)
    }

    if(length(neg_dups) > 0) {
      dup.targets[neg_dups, ]$paired.region.counts <- counts_for_ranges(dup.targets[neg_dups, ], "paired.start", "paired.end", ranges.pos.chr)
    }
    
    dup.targets$count_difference <- dup.targets[, dupcount - paired.region.counts]
    dup.targets$dups_to_keep     <- dup.targets[, ifelse(count_difference > 0, dupcount - count_difference, dupcount)]
    dup.targets$dups_to_keep     <- dup.targets[, ifelse(dups_to_keep < 1, 1, dups_to_keep)]

    dups_to_keep <- dup.targets$dups_to_keep
    dup.targets.clean <- dup.targets[, list(chr, start, end, strand)]
    kept.dups <- dup.targets.clean[rep(1:nrow(dup.targets), times=dups_to_keep), ]
    message(id, "(", chromosome, ") ", pn(nrow(kept.dups)), " duplicate reads kept out of ", pn(sum(dup.targets$dupcount)), " total")

    filtered.reads <- rbind(filtered.reads, kept.dups)
  }

  message(id, "artifact_filter: Returning ", pn(nrow(filtered.reads)), " reads after filtering")
  gr <- GRanges(seqnames=filtered.reads$chr, 
                IRanges(start=filtered.reads$start, end=filtered.reads$end), 
                strand=filtered.reads$strand,
                seqlengths=seqlengths(bam.gr))
  gr
}

counts_for_ranges <- function(df, col.start, col.end, sample.ranges) {
  df <- as.data.frame(df)
  df$counts_rowmarker <- 1:nrow(df)
  
  start.values <- df[, col.start]
  end.values   <- df[, col.end]
  
  ranged <- IRanges(start=start.values, end=end.values, names=df$counts_rowmarker)
  
  overlaps <- as.data.frame(as.matrix(findOverlaps(ranged, sample.ranges)))

  if(nrow(overlaps) == 0) {
    return(rep(0, nrow(df)))
  }
  
  ranged <- as.data.frame(ranged)
  ranged$rownum <- 1:nrow(ranged)
    
  overlaps.table <- as.data.frame(table(overlaps$queryHits))
  if(nrow(overlaps.table) > 0) {
    names(overlaps.table) <- c("rownum", "tmp_count")
    overlaps.table$rownum <- as.integer(as.character(overlaps.table$rownum))
  
    query <- merge(ranged, overlaps.table, all.x=T)
    iv.nas <- which(is.na(query$tmp_count))
    if(length(iv.nas) > 0) query$tmp_count[iv.nas] <- 0
  
    names(query)[which(names(query) == "names")] <- "counts_rowmarker"
    df <- merge(df, query[, c("counts_rowmarker", "tmp_count")])
    df <- df[order(df$counts_rowmarker), ]
    df$tmp_count
  } else {
    rep(0, nrow(df))
  }
}

pn <- function(value) {
  prettyNum(value, big.mark=",")
}


# ------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$file)) {
  message("No BAM file specified. Use --help to list available options.")
  q(status=1)
}

suppressPackageStartupMessages(library(Rsamtools, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(chipseq, warn.conflicts=F))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F))
suppressPackageStartupMessages(library(data.table, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicAlignments))

# used as a prefix for all output messages
id <- "[unknown] "

bam_file      <- opt$file
ext_length    <- opt$extension
var_name      <- opt$name
run_filter    <- !opt$skipfilter

id <- paste("[", var_name, "] ", sep="")

bigwig_file   <- opt$bigwig
paired        <- opt$paired

if(is.na(bigwig_file)) bigwig_file <- paste(var_name, "bw", sep=".")

est.frag.size <- NULL

ranges_name <- paste(var_name, ".granges", sep="")

message(id, "Converting BAM to ranges object:")
message(id, "Input BAM: ", bam_file)
message(id, "Object name: ", ranges_name)

if(!file.exists(bam_file)) {
	stop("Could not open BAM file: ", bam_file)
}

if(paired) {
  bam.gr <- readBAM_pe(bam_file)
} else {
  bam.gr <- readBAM(bam_file)
}

if(paired == FALSE & (ext_length == "auto" | run_filter == TRUE)) {
  message(id, "Estimating fragment length...")
  est.frag.size <- median(estimate.mean.fraglen(bam.gr, method="coverage"))
  message(id, "Fragment length estimate: ", est.frag.size)
}

if(paired == TRUE | ext_length == "native") {
	ext_length <- NULL
} else {
	if(ext_length == "auto")
	  ext_length <- est.frag.size
	else
	  ext_length <- as.integer(ext_length)
}

if(run_filter == TRUE & !paired) {
  bam.gr <- artifact_filter(bam.gr, ext_length=est.frag.size)
}

message(id, "Extension length: ", ifelse(is.null(ext_length), "native", ext_length))

if(!is.null(ext_length)) bam.gr <- resize(bam.gr, ext_length)

message(id, "Saving ranges object...")
bam.gr <- bam.gr[order(bam.gr)]
saveRDS(bam.gr, file=paste(ranges_name, ".rds", sep=""))
	
message(id, "Generating coverage object...")
sample.cov <- coverage(bam.gr)

nothing <- write_bigwig(sample.cov, bigwig_file)


