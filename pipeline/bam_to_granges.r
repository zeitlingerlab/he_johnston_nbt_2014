suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

option_list <- list(
  make_option(c("-f", "--file"), 
              type="character",
              help="Path of BAM file to process"))

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))

message("Reading BAM: ", opt$file)
bam.gr <- as(readGAlignments(opt$file), "GRanges")
message(" - ", pn(length(bam.gr)), " reads")

gr_file <- gsub("\\.bam", ".granges.rds", basename(opt$file))

message("Saving ranges...")                                  
saveRDS(bam.gr, file=gr_file)


