#' Opening VCF files
#'
#' @description
#' Creating an object to read VCF files line by line
#' 
#' @param filename name of file to open
#' @param regions (optional) regions to read
#'
#' @details blabla
#'
#' @returns an object of class VCFReader
#' 
#' @seealso [getFormats()]
#' @export
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' a <- openVCF(filename)
#' getFormats(a)

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

