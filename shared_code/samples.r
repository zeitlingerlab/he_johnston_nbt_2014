library(GenomicRanges)

rdata.base_path <- "samples/"

delayed_object_load_rds <- function(variable_name, filename) {
  if(!file.exists(filename)) warning("For delayed assignment of variable `", variable_name, "`: file does not exist (", filename, ")")
  delayedAssign(variable_name, 
                updateObject(readRDS(filename)), 
                assign.env = .GlobalEnv)
}

load_sample_cl <- function(sample_name, base_path=rdata.base_path) {
  var_name <- paste0(sample_name, ".cl")
  delayed_object_load_rds(var_name, paste0(base_path, sample_name, ".cl.rds"))
}

load_sample_granges <- function(sample_name, base_path=rdata.base_path) {
  var_name <- paste0(sample_name, ".gr")
  delayed_object_load_rds(var_name, paste0(base_path, sample_name, ".granges.rds"))
}

load_samples <- function() {
  samples.df <- read.csv("samples/samples.csv", stringsAsFactors=FALSE, header=TRUE)
  samples.df <- samples.df[order(samples.df$sample), ]
  
  cl_objects      <- lapply(samples.df$sample, load_sample_cl)
  granges_objects <- lapply(samples.df$sample, load_sample_granges)
  
  samples.df
}

reload_samples <- function() {
  samples.df <<- load_samples()
}

if(!exists("samples.df")) samples.df <- load_samples()

get_sample_cl <- function(sample_name) {
  get(paste0(sample_name, ".cl"))
}

get_sample_gr <- function(sample_name) {
  get(paste0(sample_name, ".gr"))
}

