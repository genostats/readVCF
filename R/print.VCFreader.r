#' @exportS3Method print VCFReader
print.VCFReader <- function(x, ...) {
  cat("A VCF Reader object with:\n")
  cat("filename: ", x$filename, "\n")
  cat("formats:", getFormats(x), "\n")
}
