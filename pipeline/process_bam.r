suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

option_list <- list(
  make_option(c("-f", "--file"), 
              type="character",
              help="Path of BAM file to process"),
  make_option(c("-n", "--name"),
              type="character",
              help="Name for resulting R object"))

opt <- parse_args(OptionParser(option_list=option_list))

message("Reading BAM: ", opt$file)
bam.gr <- as(readGAlignments(opt$file, use.names=TRUE), "GRanges")
message(" - ", pn(length(bam.gr)), " reads")

# We can't use GenomicRanges::unique here because it doesn't take into account metadata columns

message("Removing barcode duplicates...")
mcols(bam.gr)$barcode <- names(bam.gr)
names(bam.gr) <- NULL
bam.dt <- as.data.table(as.data.frame(bam.gr))
bam.dt.uniq <- unique(bam.dt)
message(" - ", pn(nrow(bam.dt.uniq)), " reads remaining")

bam.gr.uniq <- with(bam.dt.uniq, GRanges(ranges=IRanges(start=start, end=end), 
                                         seqnames=seqnames,
                                         strand=strand,
                                         seqlengths=seqlengths(bam.gr)))

message("Saving ranges...")                                  
saveRDS(bam.gr.uniq, file=paste(opt$name, ".granges.rds", sep=""))


