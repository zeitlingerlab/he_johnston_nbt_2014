
resize_around_summit <- function(gr, size) {
  trim(resize(gr, width=size, fix="center"))
}

import_and_filter_macs <- function(filename) {
  gr <- filter_chrs(import(filename))
  seqlengths(gr) <- seqlengths(Dmelanogaster)[seqlevels(gr)]
  gr <- trim(resize(gr, width=2000, fix="center"))
  gr[width(gr) == 2000]
}

distance_to_nearest_motif <- function(peaks.gr, motif.cov, sample_name) {

    motifs.gr <- as(motif.cov, "GRanges")
    motifs.gr <- subset(motifs.gr, score > 0)

    dtn <- as.data.frame(distanceToNearest(peaks.gr, motifs.gr, ignore.strand=TRUE))
    dtn.i <- data.frame(stringsAsFactors=FALSE, 
                        within_10   = length(which(dtn$distance <= 10))  / length(dtn$distance) * 100,
                        within_20   = length(which(dtn$distance <= 20))  / length(dtn$distance) * 100,
                        within_50   = length(which(dtn$distance <= 50))  / length(dtn$distance) * 100,
                        within_100  = length(which(dtn$distance <= 100)) / length(dtn$distance) * 100, 
                        sample_name = sample_name)
    dtn.i
}

distance_to_motif_plot_data <- function(chipseq, chipnexus, motif.cov) {

  random_count <- length(chipseq)
  random.gr <- sample(filter_chrs(tileGenome(tilewidth=500, seqlengths=seqlengths(Dmelanogaster), cut.last=TRUE)), random_count)
  random.gr <- resize(random.gr, width=1, fix="center")
  
  n.df <- rbind(distance_to_nearest_motif(chipseq, motif.cov, "ChIP-seq"),
                distance_to_nearest_motif(chipnexus, motif.cov, "ChIP-nexus"),
                distance_to_nearest_motif(random.gr, motif.cov, "Genome"))
    
  n.df
}

motif_profiles <- function(gr, motif.cov, window_size) {
  motif.i <- standard_metapeak(gr, motif.cov, upstream=window_size, downstream=window_size, sample_name="Motif")
  motif.i
}

reads_for_peaks <- function(gr.chipseq, gr.chipnexus, chipseq_sample, chipnexus_sample.cl, motif.cov, window_size) {

  chipseq.reads <- standard_metapeak(gr.chipseq, chipseq_sample, upstream=window_size, downstream=window_size, sample_name="ChIP-seq")
  chipseq_motifs.reads <- motif_profiles(gr.chipseq, motif.cov, window_size)
  chipseq.reads <- transform(chipseq.reads, reads = reads / max(reads) * 100)
  chipseq_motifs.reads <- transform(chipseq_motifs.reads, reads = reads * 100)

  chipnexus.reads  <- exo_metapeak(gr.chipnexus, chipnexus_sample.cl, 
                                   upstream=window_size, downstream=window_size, sample_name="ChIP-nexus", smooth=3)
  chipnexus_motifs.reads <- motif_profiles(gr.chipnexus, motif.cov, window_size)
  chipnexus.reads <- transform(chipnexus.reads, reads = reads / max(reads) * 100)
  chipnexus_motifs.reads <- transform(chipnexus_motifs.reads, reads = reads * 100)
  
  list(cs=chipseq.reads, cs_motifs=chipseq_motifs.reads,
       cn=chipnexus.reads, cn_motifs=chipnexus_motifs.reads)  
}

plots_for_factor <- function(factor_name, data.list) {
  cs.df <- data.list$cs
  cs_motifs.df <- data.list$cs_motif
  cn.df <- data.list$cn
  cn_motifs.df <- data.list$cn_motif

  cs.df <- transform(cs.df, reads = reads / max(reads) * 100)
  cn.df <- transform(cn.df, reads = reads / max(reads) * 100)

  cs.df$strand <- "ChIP-seq"
  cn.df$strand <- ifelse(cn.df$strand == "+", "+ ChIP-nexus", "- ChIP-nexus")
  profiles.df <- rbind(cs.df, cn.df)

  cs_motifs.df$sample_name <- "ChIP-seq"
  cn_motifs.df$sample_name <- "ChIP-nexus"
  motifs.df <- rbind(cs_motifs.df, cn_motifs.df)

  chipseq_color <- "#6FBE45"
  chipnexus_pos_color <- "#ED2024"
  chipnexus_neg_color <- "#2A2C7C"
  chipnexus_color <- hex(mixcolor(0.5, hex2RGB(chipnexus_pos_color), hex2RGB(chipnexus_neg_color)))

  g.profiles <- ggplot(profiles.df, aes(x=tss_distance, y=reads, color=strand)) +
                geom_line() +
                geom_area(position="identity", color=NA, aes(fill=strand), alpha=0.5) +
                geom_vline(xintercept=0, color="black", linetype="dotted") +
                scale_fill_manual("", values=c("ChIP-seq"=chipseq_color, 
                                                 "+ ChIP-nexus"=chipnexus_pos_color,
                                                 "- ChIP-nexus"=chipnexus_neg_color)) +
                scale_colour_manual("", values=c("ChIP-seq"=chipseq_color, 
                                                 "+ ChIP-nexus"=chipnexus_pos_color,
                                                 "- ChIP-nexus"=chipnexus_neg_color), guide="none") +
                theme_manuscript() +
                scale_y_continuous(expand=c(0, 0)) +
                labs(x="Distance to peak summit (bp)",
                     y="Signal (% of maximum)",
                     title=paste0(factor_name, " peak profiles"))
                   

  g.motifs <- ggplot(motifs.df, aes(x=tss_distance, y=reads, color=sample_name)) +
              geom_line() +
              geom_area(position="identity", color=NA, aes(fill=sample_name), alpha=0.5) +
              geom_vline(xintercept=0, color="black", linetype="dotted") +
              scale_fill_manual("", values=c("ChIP-seq"=chipseq_color, 
                                             "ChIP-nexus"=chipnexus_color)) +
              scale_colour_manual("", values=c("ChIP-seq"=chipseq_color, 
                                               "ChIP-nexus"=chipnexus_color)) +
              theme_manuscript() +
              scale_y_continuous(expand=c(0, 0), limits=c(0, ceiling(max(motifs.df$reads) * 1.5))) +
              labs(x="Distance to peak summit (bp)",
                   y="Motif presence (%)",
                   title=paste0(factor_name, " motif profiles"))
            
  list(profiles=g.profiles, motifs=g.motifs)
}

nearest_motif_plot <- function(factor_name, chipseq, chipnexus, motif.cov) {

  chipseq_color <- "#6FBE45"
  chipnexus_pos_color <- "#ED2024"
  chipnexus_neg_color <- "#2A2C7C"
  chipnexus_color <- hex(mixcolor(0.5, hex2RGB(chipnexus_pos_color), hex2RGB(chipnexus_neg_color)))

  near_motifs.df <- distance_to_motif_plot_data(chipseq=chipseq, chipnexus=chipnexus, motif.cov)

  motifs.m <- melt(near_motifs.df, id.var=c("sample_name"))
  motifs.m$sample_name <- factor(motifs.m$sample_name, levels=c("Genome", "ChIP-seq", "ChIP-nexus"))
  motifs.m$variable <- factor(motifs.m$variable, levels=unique(motifs.m$variable))

  g <- ggplot(motifs.m, aes(x=variable, y=value, fill=sample_name)) +
       geom_bar(stat="identity", position="dodge", color=NA) +
       scale_fill_manual(name="", values=c("ChIP-seq"=chipseq_color, 
                                           "ChIP-nexus"=chipnexus_color,
                                           "Genome"="gray50")) +
       scale_y_continuous(expand=c(0, 0)) +
       theme_manuscript() +
       labs(x="Distance to nearest motif (bp)",
            y="Percent of peaks", 
            title=factor_name)

  list(plot=g, data=motifs.m)
}

proportion_testing_results_disabled <- function(counts.df) {
  results.df <- counts.df %>% 
                  group_by(variable) %>%
                  summarize(test="ChIP-nexus > ChIP-seq", 
                            pvalue=prop.test(x=c(value[sample_name == "ChIP-nexus"], value[sample_name == "ChIP-seq"]) / 100 * 200,
                                             n=c(200, 200),
                                             alternative="greater")$p.value)
  results.df
}

chisq_testing_results <- function(counts.df) {
  results.df <- counts.df %>% 
                group_by(variable) %>%
                dplyr::summarize(test="ChIP-nexus different than ChIP-seq",
                                 pvalue=chisq.test(matrix(c(value[sample_name == "ChIP-nexus"] / 100 * 200,
                                                          200 - value[sample_name == "ChIP-nexus"] / 100 * 200,
                                                          value[sample_name == "ChIP-seq"] / 100 * 200,
                                                          200 - value[sample_name == "ChIP-seq"] / 100 * 200),
                                                          byrow=TRUE, nrow=2))$p.value)
  results.df
}



















