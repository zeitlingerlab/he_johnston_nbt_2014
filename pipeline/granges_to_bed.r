suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-f", "--file"), 
              type="character",
              default=NA,
              help="Path of GRanges object to process"),
  make_option(c("-r", "--rds"),
              action="store_true",
              default=FALSE,
              help="Input object is an RDS not RData"),
  make_option(c("-b", "--bed"),
              type="character",
              help="Output BED file"),
  make_option(c("-s", "--size"),
              type="integer",
              default=0,
              help="Resize ranges"),
  make_option(c("-i", "--index"),
              action="store_true",
              default=FALSE,
              help="Compress and index output BED"))

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$file)) {
  message("No input GRanges file specified. Use --help to list available options.")
  q(status=1)
}

suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F))

load_rdata <- function(filepath) {
  message("Loading: ", filepath)
  updateObject(get(load(filepath)))
}

load_rds <- function(filepath) {
  message("Loading: ", filepath)
  updateObject(readRDS(filepath))
}

if(opt$rds) {
  gr <- load_rds(opt$file)
} else {
  gr <- load_rdata(opt$file)
}

stopifnot(class(gr) == "GRanges")

if(opt$size != 0) {
  message("Resizing...")
  gr <- resize(gr, opt$size)
}

message("Writing: ", opt$bed)
export(gr, opt$bed, index=opt$index)

