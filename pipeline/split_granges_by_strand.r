suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-r", "--ranges"), 
              type="character",
              help="Path of GRanges file to process"))

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))

message("Loading: ", opt$ranges)
gr <- updateObject(readRDS(opt$ranges))

gr.p <- gr[strand(gr) == "+"]
gr.n <- gr[strand(gr) == "-"]

gr.p <- resize(gr.p, 1)
gr.n <- resize(gr.n, 1)

make_bigwig <- function(gr, filename, negative=FALSE) {
  gr.cov <- coverage(gr)
  if(negative) gr.cov <- gr.cov * -1
  message("Writing: ", filename)
  export(gr.cov, filename)
}

positive_file <- gsub("\\.granges\\.rds", "_positive.bw", opt$ranges)
negative_file <- gsub("\\.granges\\.rds", "_negative.bw", opt$ranges)

make_bigwig(gr.p, positive_file)
make_bigwig(gr.n, negative_file, negative=TRUE)

