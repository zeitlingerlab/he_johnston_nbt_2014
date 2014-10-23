library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)

message("Loading pentamers...")
pentamers.gr <- as(readGAlignments("pentamers.bam", use.names=TRUE), "GRanges")
message("Assigning names...")
mcols(pentamers.gr)$pentamer <- names(pentamers.gr)
names(pentamers.gr) <- NULL
message("Resizing...")
pentamers.gr <- resize(pentamers.gr, width=1, fix="center")
saveRDS(pentamers.gr, file="pentamers.gr.rds")
