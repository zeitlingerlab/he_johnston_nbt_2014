library(GenomicRanges)
library(reshape2)
library(lattice)

read_matrix <- function(gr, cov, reverse_reads=FALSE) {
  transform_function <- if(reverse_reads) { rev } else { identity }
  o <- order(gr)
  gr <- gr[o]
  rl <- as(gr, "RangesList")
  view <- RleViewsList(rleList=cov[names(rl)], rangesList=rl)
  reads.list <- viewApply(view, function(x) { transform_function(as.numeric(x)) })
  reads.m <- matrix(unlist(sapply(reads.list, as.numeric)), nrow=length(gr), byrow=TRUE)
  reads.m[o, ] <- reads.m
  reads.m
}

standard_metapeak_matrix <- function(regions.gr, sample.cov, upstream=100, downstream=100) {
  regions.gr <- resize(regions.gr, width=downstream)
  regions.gr <- resize(regions.gr, width=upstream + width(regions.gr), fix="end")
  
  stopifnot(length(unique(width(regions.gr))) == 1)
  
  reads <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1])
  
  i_p <- which(strand(regions.gr) == "+" | strand(regions.gr) == "*")
  i_n <- which(strand(regions.gr) == "-")
  
  if(class(sample.cov) == "character") {
    sample.cov <- import.bw(sample.cov, which=regions.gr, as="RleList")
  }
  
  if(length(i_p) > 0) reads[i_p, ] <- read_matrix(regions.gr[i_p], sample.cov)
  if(length(i_n) > 0) reads[i_n, ] <- read_matrix(regions.gr[i_n], sample.cov, reverse_reads=TRUE)

  reads
}

exo_metapeak_matrix <- function(regions.gr, sample.cl, upstream=100, downstream=100) {
  regions.gr <- resize(regions.gr, width=downstream)
  regions.gr <- resize(regions.gr, width=upstream + width(regions.gr), fix="end")

  i_p <- which(strand(regions.gr) == "+" | strand(regions.gr) == "*")
  i_n <- which(strand(regions.gr) == "-")
  
  reads.p <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1])
  reads.n <- reads.p
  
  
  if(length(i_p) > 0) {
    reads.p[i_p, ] <- read_matrix(regions.gr[i_p], sample.cl$pos)
    reads.n[i_p, ] <- read_matrix(regions.gr[i_p], sample.cl$neg)
  }
  
  if(length(i_n) > 0) {
    reads.p[i_n, ] <- read_matrix(regions.gr[i_n], sample.cl$neg, reverse_reads=TRUE)
    reads.n[i_n, ] <- read_matrix(regions.gr[i_n], sample.cl$pos, reverse_reads=TRUE)
  }

  list(pos=reads.p, neg=reads.n)
}

standard_metapeak <- function(gr, sample.cov, upstream=100, downstream=100, sample_name=NA, smooth=NA) {
  message("standard metapeak: ", sample_name)
  
  reads <- standard_metapeak_matrix(gr, sample.cov, upstream, downstream)
  
  reads.df <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                        reads=colMeans(reads), 
                        sample_name=sample_name)
  if(!is.na(smooth)) reads.df$reads <- as.numeric(runmean(Rle(reads.df$reads), k=smooth, endrule="constant"))
  reads.df  
}

exo_metapeak <- function(gr, sample.cl, upstream=100, downstream=100, sample_name=NA, smooth=NA) {
  message("exo metapeak: ", sample_name)
  reads.list <- exo_metapeak_matrix(gr, sample.cl, upstream, downstream)
  
  reads.p <- reads.list$pos
  reads.n <- reads.list$neg
  
  df.p <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                     reads=colMeans(reads.p), 
                     strand="+")

  df.n <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                     reads=colMeans(reads.n), 
                     strand="-")

  if(!is.na(smooth)) {
    df.n$reads <- as.numeric(runmean(Rle(df.n$reads), k=smooth, endrule="constant"))
    df.p$reads <- as.numeric(runmean(Rle(df.p$reads), k=smooth, endrule="constant"))
  }

  reads.df <- rbind(df.p, df.n)
  reads.df$sample_name <- sample_name
  reads.df  
}

base_frequencies <- function(gr, upstream, downstream, genome=Dmelanogaster) {

  gr <- resize(gr, width=downstream)
  gr <- resize(gr, width=upstream + width(gr), fix="end")

  m <- consensusMatrix(getSeq(genome, gr))[1:4, ]
  m <- m / colSums(m)
  df.bases <- as.data.frame(t(m))
  df.bases$tss_distance <- (-1 * upstream):(downstream - 1)
  df.bases <- melt(df.bases, id.var="tss_distance")
  names(df.bases)[2:3] <- c("base", "frequency")
  df.bases
}

normalize_matrix <- function(m, value) {
  m <- pmin(m / value, 1)
  m
}

draw_standard_heatmap <- function(m, title, center_vline=TRUE, normalize=TRUE) {
  pos.colors <- colorRampPalette(c("white", "red"))(32)
  
  if(normalize) {
    max.value <- quantile(m, 0.98)
    m <- normalize_matrix(m, max.value)
  }

  image(t(m[nrow(m):1, ]), col=pos.colors, useRaster=TRUE, main=title, yaxt="n", xaxt="n")

  if(center_vline) {
    vline <- matrix(NA, nrow=nrow(m), ncol=ncol(m))
    
    m_center_col <- ceiling(ncol(m) / 2)
    m_line_width <- ceiling(ncol(m) * 0.01)
    m_center_range <- (m_center_col - m_line_width) : (m_center_col + m_line_width)
    vline[, m_center_range] <- matrix(1, nrow=nrow(vline), ncol=length(m_center_range))
    vline.colors <- paste0(colorRampPalette("gray")(1), "77")
    image(t(vline[nrow(vline):1,]), col=vline.colors, useRaster=TRUE, add=TRUE, yaxt="n", xaxt="n")
  }
}

draw_exo_heatmap <- function(strand.list, title, sort_column_range=1:ncol(m.pos)) {
  pos.colors <- paste0(colorRampPalette(c("white", "#D62A31"))(33), "88")
  neg.colors <- paste0(colorRampPalette(c("white", "#1E214A"))(33), "88")
  
  m.pos <- strand.list$pos
  m.neg <- strand.list$neg
  
  stopifnot(ncol(m.pos) == ncol(m.neg))
  stopifnot(nrow(m.pos) == nrow(m.neg))
  
  max.value <- quantile(rbind(m.pos, m.neg), 0.98)
  
  m.pos <- normalize_matrix(m.pos, max.value)
  m.neg <- normalize_matrix(m.neg, max.value)

  m.both <- m.pos + m.neg
  
  row.order <- order(rowSums(m.both[, sort_column_range]), decreasing=TRUE)
  
  m.pos <- m.pos[row.order, ]
  m.neg <- m.neg[row.order, ]
  
  m.pos[m.pos == 0] <- NA
  m.neg[m.neg == 0] <- NA
  
  # Combined image
  image(t(m.pos[nrow(m.pos):1,]), col=pos.colors, useRaster=TRUE, main=title)
  image(t(m.neg[nrow(m.neg):1,]), col=neg.colors, useRaster=TRUE, add=TRUE)
  
  # Separate image (for scales)
  print(levelplot(t(m.pos[nrow(m.pos):1, ]), col.regions=pos.colors, useRaster=TRUE, main=paste0(title, " - positive strand"), ylab="", xlab="", cuts=31))
  print(levelplot(t(m.neg[nrow(m.neg):1, ]), col.regions=neg.colors, useRaster=TRUE, main=paste0(title, " - negative strand"), ylab="", xlab="", cuts=31))
}

