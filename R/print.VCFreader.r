#' Print a VCFReader object
#'
#' @description
#' Prints a short summary of a `VCFReader` object, including the filename and available FORMAT fields.
#'
#' @param x An object of class `VCFReader`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return No return value. This function is invoked for its side effect (printing to the console).
#'
#' @seealso \link{VCFReader}, \link{openVCF}
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' reader <- openVCF(filename)
#' print(reader)  # or just `reader`
#'
#'
#' @exportS3Method print VCFReader
print.VCFReader <- function(x, ...) {
  cat("A VCF Reader object with:\n")
  cat("filename: ", x$filename, "\n")
  cat("formats:", getFormats(x), "\n")
}
