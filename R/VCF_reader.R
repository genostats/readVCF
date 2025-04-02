openVCF <- function(filename, regions = c()) {
  ptr <- .Call("_readVCF_openVCF", filename, regions) 
  structure(ptr, class = "VCFReader")
}

print <- function(x, ...) {
  cat("A VCF Reader\n")
}