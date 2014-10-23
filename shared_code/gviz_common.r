
gviz_plot_with_exo <- function(genome="dm3", region.gr, 
                               sample.regular, sample.nexus, sample.exo, 
                               title, motifs=GRanges(),
                               normalization_method="none",
                               zoom_to_nexus_scale=FALSE) {

  r_chr <- as.character(seqnames(region.gr))
  r_start <- start(region.gr)
  r_end   <- end(region.gr)

  stopifnot(length(region.gr) == 1)

  seqlevels(region.gr, force=TRUE) <- unique(as.character(seqnames(region.gr)))

  if(genome == "dm3") {
    grtrack <- GeneRegionTrack(range=TxDb.Dmelanogaster.UCSC.dm3.ensGene,
                               genome="dm3",
                               chromosome=r_chr,
                               name="Genes",
                               showId=TRUE,
                               geneSymbol=TRUE,
                               fill="gray60", collapseTranscripts=FALSE) 
  }
  
  if(genome == "hg19") {
    grtrack <- GeneRegionTrack(range=TxDb.Hsapiens.UCSC.hg19.knownGene,
                               genome="hg19",
                               chromosome=r_chr,
                               name="Genes",
                               showId=TRUE,
                               geneSymbol=TRUE,
                               fill="gray60", collapseTranscripts=FALSE) 
  }

  if(normalization_method == "downsample") {
    if(length(sample.nexus$gr) > length(sample.exo$gr)) {
      sample.nexus$gr <- sample(sample.nexus$gr, length(sample.exo$gr))
      sample.nexus$pos <- coverage(sample.nexus$gr[strand(sample.nexus$gr) == "+"])
      sample.nexus$neg <- coverage(sample.nexus$gr[strand(sample.nexus$gr) == "-"])
    } else {
      sample.exo$gr <- sample(sample.exo$gr, length(sample.nexus$gr))
      sample.exo$pos <- coverage(sample.exo$gr[strand(sample.exo$gr) == "+"])
      sample.exo$neg <- coverage(sample.exo$gr[strand(sample.exo$gr) == "-"])
    }
  }

  

  data.nexus <- matrix(c(as.numeric(sample.nexus$pos[[r_chr]][r_start:r_end]),
                         -1 * as.numeric(sample.nexus$neg[[r_chr]][r_start:r_end])), nrow=2, byrow=TRUE)
  
  data.exo <- matrix(c(as.numeric(sample.exo$pos[[r_chr]][r_start:r_end]),
                     -1 * as.numeric(sample.exo$neg[[r_chr]][r_start:r_end])), nrow=2, byrow=TRUE)

  if(normalization_method == "read_count") {
    nexus.rc <- length(sample.nexus$gr)
    exo.rc   <- length(sample.exo$gr)
    target.rc <- ceiling(mean(c(nexus.rc, exo.rc)))
    
    data.nexus <- data.nexus / nexus.rc * target.rc
    data.exo   <- data.exo / exo.rc * target.rc
  }

  if(normalization_method == "downscale_exo") {
    nexus.rc <- length(sample.nexus$gr)
    exo.rc   <- length(sample.exo$gr)
    
    data.exo   <- data.exo / exo.rc * nexus.rc
  }


  ylim.nexus <- max(abs(c(max(data.nexus), min(data.nexus))))
  ylim.exo <- max(abs(c(max(data.exo), min(data.exo))))

  if(zoom_to_nexus_scale) ylim.exo <- ylim.nexus

  nexus_track <- DataTrack(data=data.nexus,
                         start=r_start:r_end, width=0, chromosome=r_chr,
                         genome=genome, name="ChIP-nexus", 
                         groups=c("Positive strand", "Negative strand"), 
                         type="histogram",
                         ylim=c(-ylim.nexus, ylim.nexus),
                         legend=TRUE,
                         col=c("darkblue", "red"), col.line=c("darkblue", "red"))

  exo_track <- DataTrack(data=data.exo,
                         start=r_start:r_end, width=0, chromosome=r_chr,
                         genome=genome, name="ChIP-exo", 
                         groups=c("Positive strand", "Negative strand"), 
                         type="histogram",
                         ylim=c(-ylim.exo, ylim.exo),
                         legend=TRUE,
                         col=c("darkblue", "red"), col.line=c("darkblue", "red"))

  regular.rle <- import(sample.regular, which=region.gr, as="RleList")
  data.reg <- matrix(c(as.numeric(regular.rle[[r_chr]][r_start:r_end])), nrow=1, byrow=TRUE) 
  ylim.reg <- max(data.reg)
  
  reg_track <- DataTrack(range=sample.regular,
                             genome=genome, name="ChIP-seq", 
                             type="polygon",
                             ylim=c(0, ylim.reg),
                             col=c("#1C5425"),
                             fill.mountain=c("#1C5425", "#1C5425"))


  gtrack <- GenomeAxisTrack()

  tlist <- list(gtrack, reg_track, nexus_track, exo_track, grtrack)      
  tsizes <- c(0.1, 0.25, 1, 1, 0.1)

  if(length(motifs) > 0) {
    motifs <- motifs[seqnames(motifs) == r_chr]
    motif_track <- AnnotationTrack(range=motifs, strand=rep("*", length(motifs)),
                                  genome=genome, name="Motifs", showFeatureId=TRUE, 
                                  stacking="dense", fill="gray50", fontcolor="black", fontcolor.group="black", fontsize=5)
    tlist <- c(tlist, list(motif_track))
    tsizes <- c(tsizes, 0.1)
  }

  plotTracks(tlist, 
             sizes=tsizes,
             chromosome=r_chr, 
             from=r_start,
             to=r_end,
             main=title,
             cex.title=1.2, cex.axis=0.8, col.title="black", col.axis="black",
             fontcolor.legend="black", cex.legend=1.1) 
}


gviz_plot_no_exo <- function(genome="dm3", region.gr, 
                             sample.regular, sample.nexus, 
                             title, motifs=GRanges()) {

  r_chr <- as.character(seqnames(region.gr))
  r_start <- start(region.gr)
  r_end   <- end(region.gr)

  stopifnot(length(region.gr) == 1)

  seqlevels(region.gr, force=TRUE) <- unique(as.character(seqnames(region.gr)))

  if(genome == "dm3") {
    grtrack <- GeneRegionTrack(range=TxDb.Dmelanogaster.UCSC.dm3.ensGene,
                               genome="dm3",
                               chromosome=r_chr,
                               name="Genes",
                               showId=TRUE,
                               geneSymbol=TRUE,
                               fill="gray60", collapseTranscripts=FALSE) 
  }
  
  if(genome == "hg19") {
    grtrack <- GeneRegionTrack(range=TxDb.Hsapiens.UCSC.hg19.knownGene,
                               genome="hg19",
                               chromosome=r_chr,
                               name="Genes",
                               showId=TRUE,
                               geneSymbol=TRUE,
                               fill="gray60", collapseTranscripts=FALSE) 
  }

  data.nexus <- matrix(c(as.numeric(sample.nexus$pos[[r_chr]][r_start:r_end]),
                         -1 * as.numeric(sample.nexus$neg[[r_chr]][r_start:r_end])), nrow=2, byrow=TRUE)
  

  ylim.nexus <- max(abs(c(max(data.nexus), min(data.nexus))))

  nexus_track <- DataTrack(data=data.nexus,
                         start=r_start:r_end, width=0, chromosome=r_chr,
                         genome=genome, name="ChIP-nexus", 
                         groups=c("Positive strand", "Negative strand"), 
                         type="histogram",
                         ylim=c(-ylim.nexus, ylim.nexus),
                         legend=TRUE,
                         col=c("darkblue", "red"), col.line=c("darkblue", "red"))

  regular.rle <- import(sample.regular, which=region.gr, as="RleList")
  data.reg <- matrix(c(as.numeric(regular.rle[[r_chr]][r_start:r_end])), nrow=1, byrow=TRUE) 
  ylim.reg <- max(data.reg)
  
  reg_track <- DataTrack(range=sample.regular,
                             genome=genome, name="ChIP-seq", 
                             type="polygon",
                             ylim=c(0, ylim.reg),
                             col=c("#1C5425"),
                             fill.mountain=c("#1C5425", "#1C5425"))


  gtrack <- GenomeAxisTrack()

  tlist <- list(gtrack, reg_track, nexus_track, grtrack)      
  tsizes <- c(0.1, 0.25, 1, 0.1)

  if(length(motifs) > 0) {
    motifs <- motifs[seqnames(motifs) == r_chr]
    motif_track <- AnnotationTrack(range=motifs, strand=rep("*", length(motifs)),
                                  genome=genome, name="Motifs", showFeatureId=TRUE, 
                                  stacking="dense", fill="gray50", fontcolor="black", fontcolor.group="black", fontsize=5)
    tlist <- c(tlist, list(motif_track))
    tsizes <- c(tsizes, 0.1)
  }

  plotTracks(tlist, 
             sizes=tsizes,
             chromosome=r_chr, 
             from=r_start,
             to=r_end,
             main=title,
             cex.title=1.2, cex.axis=0.8, col.title="black", col.axis="black",
             fontcolor.legend="black", cex.legend=1.1) 
}



