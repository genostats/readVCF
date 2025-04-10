openVCF <- function(filename, regions = character(0)) {
  if(length(regions) == 0) ptr <- .Call("_readVCF_openVCFonly", filename) 
  else ptr <- .Call("_readVCF_openVCFregs", filename, regions) 
  structure(ptr, class = "VCFReader")
}

#x = 'VCFReader' in lieu of objects ?
print.VCFReader <- function(obj='VCFReader') {
  cat("A VCF Reader\n")
  #ptr->regions
}
