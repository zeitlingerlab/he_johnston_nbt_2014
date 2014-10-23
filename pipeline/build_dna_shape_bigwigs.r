library(GenomicRanges)
library(rtracklayer)

pentamers.gr <- readRDS("pentamers.gr.rds")

shape.df <- readRDS("shape.df.rds")

build_bigwig <- function(pentamers.gr, shape.df, parameter, bigwig_file) {
  nothing <- gc()
  message(parameter)
  message(" - assigning values")
  mcols(pentamers.gr)$shape_param <- shape.df[, parameter][match(mcols(pentamers.gr)$pentamer, shape.df$pentamer)]
  message(" - calculating coverage")
  shape.cov <- coverage(pentamers.gr, weight="shape_param")
  message(" - writing bigwig")
  export(shape.cov, bigwig_file)
}

build_bigwig(pentamers.gr, shape.df, "mgw", "minor_groove_width.bw")
build_bigwig(pentamers.gr, shape.df, "prop_twist", "propeller_twist.bw")


