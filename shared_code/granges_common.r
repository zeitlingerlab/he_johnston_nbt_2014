library(GenomicRanges)
library(rtracklayer)

check_coverage_argument <- function(cvg, regions=NULL) {
  if(class(cvg) == "character") {
    if(is.null(regions)) {
      cvg <- import(cvg, as="RleList")
    } else {
      cvg <- import(cvg, which=regions, as="RleList")
    }
  } 
  cvg
}

regionSums <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewSums(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMeans <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMeans(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionWhichMaxs <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewWhichMaxs(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMaxs <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMaxs(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMins <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMins(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

total_signal <- function(cov) {
  cov <- check_coverage_argument(cov)
  sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x)))))
}

get_rds <- function(filename) {
  message("Loading ", filename, " ... ", appendLF=FALSE)
  o <- updateObject(readRDS(filename))
  message("OK")
  o
}

get_load <- function(filename) {
  message("Loading ", filename, " ... ", appendLF=FALSE)
  o <- updateObject(get(load(filename)))
  message("OK")
  o
}

filter_chrs <- function(gr) {
  keep.chrs <- seqlevels(gr)[-grep("H|U|M|_", seqlevels(gr))]
  if(length(keep.chrs) > 0) seqlevels(gr, force=TRUE) <- keep.chrs
  gr
}
