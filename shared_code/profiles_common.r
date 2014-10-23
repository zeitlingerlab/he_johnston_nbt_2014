
top_motifs_for_sample <- function(sample_name, motifs.gr, n=200, region_size=15) {
  motif_regions.gr <- resize(motifs.gr, width=region_size, fix="center")
  mcols(motifs.gr)$signal <- regionSums(motif_regions.gr, get_sample_cl(sample_name)$cov)
  
      
  motifs.no_strand <- motifs.gr
  strand(motifs.no_strand) <- "*"

  motif.areas <- as(slice(coverage(reduce(motifs.no_strand)), 1), "GRanges")

  ol <- as.data.frame(findOverlaps(motif.areas, motifs.no_strand, select="all", ignore.strand=TRUE))
  ol$exo_signal <- mcols(motifs.gr)$signal[ol$subjectHits]
  
  ol <- arrange(ol, queryHits, plyr::desc(exo_signal))
  ol <- ol[!duplicated(ol$queryHits), ]

  stopifnot(length(which(duplicated(ol$subjectHits) == FALSE)) == nrow(ol))

  motifs.gr <- motifs.gr[ol$subjectHits]
  
  motifs.gr <- motifs.gr[order(mcols(motifs.gr)$signal, decreasing=TRUE)]
  motifs.gr[1:n]
}

process_sample <- function(sample_name, motifs.gr, n=200, region_size=15) {
  top_motifs.gr <- top_motifs_for_sample(sample_name, motifs.gr, n=n, region_size=region_size)
  export(top_motifs.gr, figure_path(paste0(sample_name, "_top_", length(top_motifs.gr), "_motifs.bed")))

  top_motifs.gr <- resize(top_motifs.gr, width=1)

  reads <- exo_metapeak(top_motifs.gr, get_sample_cl(sample_name), upstream=40, downstream=50)
  
  reads$sample_name <- sample_name
  reads$motif_width <- width(motifs.gr)[1]
  reads
}

build_plot <- function(reads.df) {

  plot_title <- paste0(reads.df$sample_name[1], " (", subset(samples.df, sample == reads.df$sample_name[1])$digestion, ")")

  motif.box <- data.frame(xmin=0, 
                          xmax=reads.df$motif_width[1]-1,
                          ymin=-Inf,
                          ymax=Inf)
                          
  g <- ggplot(reads.df, aes(x=tss_distance, y=reads, color=strand)) +
       geom_line() +
       geom_rect(show_guide=FALSE, inherit.aes=FALSE, data=motif.box, 
                 aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
                 alpha=0.25, fill="gray80") +
       theme_bw() +
       labs(x="Distance to motif left edge", y="Average reads", title=plot_title)
  g
}

heatmap_reads <- function(bedfile) {
  sample_name <- gsub("_top_200_motifs.bed", "", bedfile)
  cache.file <- paste0(sample_name, ".heatmap.matrix.rds")
  cache(cache.file, function() {
    gr <- import(figure_path(bedfile), asRangedData=FALSE)
    exo_metapeak_matrix(gr, get_sample_cl(sample_name), upstream=100, downstream=100 + 9)
  })
}