openVCF <- function(filename, field = c("GT", "DS"), regions = character(0)) {
  field <- match.arg(field)
  if(field == "GT")
    ptr <- .Call("_readVCF_openVCF_GT", filename, regions) 
  else if(field == "DS")
    ptr <- .Call("_readVCF_openVCF_DS", filename, regions) 

  structure(ptr, class = "VCFReader")
}

#x = 'VCFReader' in lieu of objects ?
print.VCFReader <- function(obj='VCFReader') {
  cat("A VCF Reader\n")
  #ptr->regions
}
