library(digest)

cache <- function(filename, func) {
  filename <- figure_path(filename)
  hashfile <- paste0(filename, ".hash")
  current_hash <- digest(paste0(as.character(body(func)), collapse="\n"), NULL, ascii=TRUE)
  new_hash_created <- FALSE
  
  if(!file.exists(hashfile)) {
    new_hash_created <- TRUE
    saveRDS(current_hash, file=hashfile)
  }
  
  if(file.exists(filename)) {
    saved_hash <- readRDS(hashfile)
    if(saved_hash != current_hash) {
      warning("For cache file '", filename, "': function has been modified since results were cached.")
      cache_results <- func()
      saveRDS(cache_results, file=filename)
      saveRDS(current_hash, file=hashfile)
      return(cache_results)
    }
    if(new_hash_created) {
      warning("For cache file '", filename, "': function hash was not present.")
    }
    return(readRDS(filename))
  } else {
    cache_results <- func()
    saveRDS(cache_results, file=filename)
    saveRDS(current_hash, file=hashfile)
    return(cache_results)
  }
}

clear_caches <- function() {
  file.remove(list.files(figure_path(), full.names=TRUE))
}
