
checked_mclapply <- function(...) {
  results <- mclapply(...)
  errors <- which(sapply(results, class) == "try-error")
  if(length(errors) > 0) stop("mclapply() returned errors: ", paste0(results[errors], collapse=", "))
  results
}
