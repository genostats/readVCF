#' Compare dosage matrices from DS and GP fields
#'
#' @description
#' This function compares the dosage matrices obtained from a VCF file using the `"DS"` and `"GP"` FORMAT fields.
#' It generates the corresponding `.dosf` files using `mk.dosage.file()`, then loads and compares them.
#'
#' @param vcf.file Path to a VCF file (preferably bgzipped and indexed).
#' @param regions Optional character vector of genomic regions to extract (e.g., `"19:320000-330000"`). Defaults to reading the entire file.
#' @param tol Numeric tolerance threshold. If the maximum absolute difference between the two matrices exceeds this value, an error is raised.
#'
#' @details
#' \itemize{
#'   \item Extracts dosage values from both `"DS"` and `"GP"` fields using `mk.dosage.file()`.
#'   \item Loads the resulting `.dosf` matrices with [bigmemory::attach.big.matrix()].
#'   \item Compares the dimensions and values of the matrices.
#'   \item If they are equal within the given tolerance, prints a success message and previews a 5x5 submatrix.
#'   \item Otherwise, an error is raised.
#' }
#'
#' @return Returns `invisible(TRUE)` if the matrices match within the tolerance.
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' test.dosage.field.match(filename, regions = "19:320000-330000")
#'
#' @seealso \link{mk.dosage.file}, \link[bigmemory]{attach.big.matrix}
test.dosage.field.match <- function(vcf.file, regions = NULL, tol = 1e-3) {
  stopifnot(file.exists(vcf.file))
  if (is.null(regions)) regions <- character(0)
  
  # Generation of files with field = 'DS' and 'GP'
  prefix_ds <- tempfile("dosage_ds_")
  prefix_gp <- tempfile("dosage_gp_")
  
  mk.dosage.file(vcf.file, prefix_ds, regions = regions, type = "float", field = "DS")
  mk.dosage.file(vcf.file, prefix_gp, regions = regions, type = "float", field = "GP")
  
  # Loading the matrices
  library(bigmemory)
  dos_ds <- attach.big.matrix(paste0(prefix_ds, ".dosf.desc"))
  dos_gp <- attach.big.matrix(paste0(prefix_gp, ".dosf.desc"))
  
  # Comparison of matrix dimensions
  if (!all(dim(dos_ds) == dim(dos_gp))) {
    stop("Matrices have different dimensions: DS = ", paste(dim(dos_ds), collapse = "x"),
         ", GP = ", paste(dim(dos_gp), collapse = "x"))
  }
  
  cat("Aperçu des premières valeurs (max 5x5) :")
  rmax <- min(5, nrow(dos_ds))
  cmax <- min(5, ncol(dos_ds))
  cat("\nDS matrix:")
  print(dos_ds[1:rmax, 1:cmax])
  cat("GP matrix:")
  print(dos_gp[1:rmax, 1:cmax])
  
  cat(paste0("\nComparaison des valeurs avec tolérance de ", tol))
  diff <- abs(dos_ds[,] - dos_gp[,])
  max_diff <- max(diff, na.rm = TRUE)
  
  if (max_diff > tol) {
    stop("\nLes matrices diffèrent au-delà de la tolérance : max diff = ", max_diff)
  } else {
    cat(paste0("\nLes matrices sont équivalentes à une tolérance près (max diff = ", max_diff, ").\n"))
  }
    
  invisible(TRUE)
}      
