suppressPackageStartupMessages(library(GenomicRanges))

files <- list.files(".", "granges.rds")

for(f in files) {
  message(f, " ", appendLF=FALSE)
  f.gr <- updateObject(readRDS(f))
  message(length(f.gr))
}

