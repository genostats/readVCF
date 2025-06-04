#' 
#' @param filename
#' @param regions
#'
#' @export 
openVCF <- function(filename, regions) {
  if(missing(regions))
    regions <- character(0)

  fi <- path.expand(filename)
  if(!file.exists(fi))
    stop("File doesn't exist")

  L <- list(filename = filename, xptr = openVCFregions(fi, regions))
  class(L) <- "VCFReader"
  L
}

