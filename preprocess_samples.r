library(GenomicRanges)
library(parallel)
library(yaml)

output_path <- "samples/"
original_rdata_path <- file.path(yaml.load_file("_output.yaml")$data$base_path, "rdata/")

samples.df <- read.csv("samples/samples.csv", stringsAsFactors=FALSE, header=TRUE)

process_sample <- function(i, df) {
  message(df$sample[i])
  ranges.file <- paste0(output_path, df$sample[i], ".granges.rds")
  if(file.exists(ranges.file)) {
    message(" - already processed")
  } else {
    message(" - Loading original ranges")
    gr <- updateObject(readRDS(paste0(original_rdata_path, df$sample[i], ".granges.rds")))

    message(" - Ordering ranges")
    gr <- gr[order(gr)]

    message(" - Resizing")
    gr <- resize(gr, width=1)
    message(" - Saving ", ranges.file)
    saveRDS(gr, file=ranges.file)

    message(" - Building total coverage")
    gr.cov <- coverage(gr)

    message(" - Building + strand coverage")
    pos.cov <- coverage(gr[strand(gr) == "+"])

    message(" - Building - strand coverage")
    neg.cov <- coverage(gr[strand(gr) == "-"])
    
    cl <- list(gr=gr, cov=gr.cov, pos=pos.cov, neg=neg.cov)

    cl.file <- paste0(output_path, df$sample[i], ".cl.rds")

    message(" - Saving ", cl.file)
    saveRDS(cl, file=cl.file)
  }
  TRUE
}

nothing <- mclapply(1:nrow(samples.df), process_sample, samples.df, mc.cores=4)



