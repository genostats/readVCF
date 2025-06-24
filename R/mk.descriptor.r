#' Create a descriptor file for a file-backed big.matrix
#'
#' @description
#' Generates a `.desc` descriptor file compatible with `bigmemory::attach.big.matrix()`,
#' enabling future access to a binary file-backed matrix (such as a `.dosf` file).
#'
#' This is particularly useful when storing large matrices on disk and reattaching them
#' later without reloading the data in memory.
#'
#' @param path A character string giving the base path (without extension) for the descriptor.
#' It is typically the path to the `.dosf` file (e.g., `"/tmp/matrix.dosf"`).
#' @param nrow An integer: number of rows in the matrix.
#' @param ncol An integer: number of columns in the matrix.
#' @param type A character string specifying the storage type: one of `"char"`, `"short"`, `"integer"`, `"double"`, or `"float"`.
#'
#' @details
#' The function creates a descriptor `.desc` file corresponding to a file-backed matrix.
#' It assumes that the associated binary data file (`.dosf`) already exists.
#' The matrix is assumed to be stored in **column-major order** (i.e., filled column by column).
#'
#' If the descriptor file already exists, it will not be overwritten, and a warning is issued.
#' To open the matrix later in R, use [bigmemory::attach.big.matrix()].
#'
#' This function is especially useful after generating genotype dosage matrices from VCF files using [mk.dosage.file()].
#'
#' @return `invisible(NULL)`. Called for its side effect (writing a `.desc` file).
#'
#' @seealso \link[bigmemory]{attach.big.matrix}, \link[bigmemory]{filebacked.big.matrix}, \link{mk.dosage.file}
#'
#' @examples
#' 
#' # Simulate binary matrix file with float values
#' path <- tempfile()
#' binfile <- paste0(path, ".dosf")
#'
#' # Write a 4x5 matrix (column-major) to .dosf
#' con <- file(binfile, "wb")
#' writeBin(as.numeric(1:20), con, size = 4)  # 4 bytes for float
#' close(con)
#'
#' # Create the corresponding descriptor
#' mk.descriptor.file(path = binfile, nrow = 4, ncol = 5, type = "float")
#'
#' # Reopen the matrix and display it
#' dosages <- bigmemory::attach.big.matrix(paste0(binfile, ".desc"))
#' print(dosages[, ])
#' 
#'
#' @export
mk.descriptor.file <- function(path, nrow, ncol, type) {
  dir <- dirname(path)
  fil <- basename(path)
  
  d <- sprintf("new(\"big.matrix.descriptor\", description = list(sharedType = \"FileBacked\",\n filename = \"%s\", ", fil)
  d <- paste0(d, sprintf("dirname = \"%s\",\n ", dir))
  d <- paste0(d, sprintf("totalRows = %dL, totalCols = %dL,\n ", nrow, ncol))
  d <- paste0(d, sprintf("rowOffset = c(0, %dL), colOffset = c(0, %dL),\n ", nrow, ncol))
  d <- paste0(d, sprintf("nrow = %d, ncol = %d,\n ", nrow, ncol))
  d <- paste0(d, sprintf("rowNames = NULL, colNames = NULL, type = \"%s\", separated = FALSE))\n", type))

  desc.file <- paste0(path, ".desc")
  if(file.exists(desc.file)) {
    warning(desc.file, " already exists, won't erase")
    return(invisible(NULL))
  }
  cat(d, file = desc.file)
}


