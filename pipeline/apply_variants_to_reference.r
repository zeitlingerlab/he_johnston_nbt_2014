suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

option_list <- list(
  make_option(c("-f", "--file"), 
              type="character",
              help="Path of VCF file to process"))

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm3))
suppressPackageStartupMessages(library(VariantAnnotation))

message("Loading: ", opt$file)
vcf.gr <- rowData(readVcf(opt$file, genome="dm3"))
vcf.gr <- vcf.gr[elementLengths(mcols(vcf.gr)$ALT) == 1]

message(" ", pn(length(vcf.gr)), " total variants to apply")

dm3.gr <- GRanges(ranges   = IRanges(start=1, end=seqlengths(Dmelanogaster)),
                  seqnames = names(seqlengths(Dmelanogaster)),
                  strand   = "+")
dm3.seq <- getSeq(Dmelanogaster, dm3.gr)
names(dm3.seq) <- seqnames(dm3.gr)

new_reference <- list()

for(i in seq_along(dm3.seq)) {
  chr.name <- names(dm3.seq)[i]
  message(" ", chr.name)
  chr.seq <- dm3.seq[[i]]
  chr.gr <- vcf.gr[seqnames(vcf.gr) == chr.name]
  if(length(chr.gr) == 0) {
    message("  - no variants in this chromosome")
    next
  }
  message("  - ", length(chr.gr), " variants in this chromosome")

  ranges.gr <- IRanges(start=start(chr.gr), end=end(chr.gr))
  alt.dna   <- as.character(unlist(mcols(chr.gr)$ALT))

  chr.new <- replaceAt(chr.seq, at=ranges.gr, value=alt.dna)
  new_reference <- c(new_reference, list(chr.new))
}

new_reference <- DNAStringSet(new_reference)
names(new_reference) <- names(dm3.seq)

writeXStringSet(new_reference, file="dm3_oregonr.fasta")
