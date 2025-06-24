#' Open a VCF file for line-by-line reading
#'
#' @description
#' Opens a VCF file and returns a `VCFReader` object that allows efficient, low-level access to its contents line by line.
#' This is particularly useful for large VCF files where reading everything into memory is not feasible.
#'
#' @param filename A character string specifying the path to the VCF file (can be `.vcf` or `.vcf.gz`).
#'                 The file must be tabix-indexed if regions are specified.
#' @param regions An optional character vector of genomic regions to extract (e.g., `"19:320000-700000"`). 
#'                If missing, the entire file is accessible.
#'
#' @details
#' The returned object is of class `VCFReader` and can be used with various functions from the package to access VCF content.
#'
#' @return An object of class `VCFReader` containing:
#' \describe{
#'   \item{filename}{The original file path.}
#'   \item{xptr}{An external pointer to the internal C++ reader.}
#' }
#'
#' @seealso \link{VCFReader}
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' reader <- openVCF(filename)
#' getFormats(reader)  # e.g., "GT", "DS", "AD", ...
#'
#' # Iterate over a region
#' reader <- openVCF(filename, "19:320000-700000")
#' while (VCFnext(reader)) {
#'   line <- getLine(reader, "DS")
#'   print(head(line$values))
#' }
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

