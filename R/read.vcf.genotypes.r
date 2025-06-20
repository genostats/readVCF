#' Read genotype matrix from a VCF file
#'
#' @description
#' Reads genotype data from a VCF file and returns it as an integer matrix. By default, all variants in the file are read,
#' unless a subset of genomic regions is specified..
#' @param filename A character string giving the path to the VCF file (can be gzipped and tabix-indexed).
#' @param regions An optional character vector of genomic regions to extract (e.g., `"19:320000-700000"`).
#' If omitted, the entire file is processed.
#'
#' @return An integer matrix with:
#' \describe{
#'   \item{Rows:}{Sample IDs (individuals).}
#'   \item{Columns:}{Variant IDs (e.g., rsIDs from the VCF).}
#'   \item{Values:}{Genotypes coded as integers (typically 0, 1, 2). Missing values may appear as `NA`.}
#' }
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' genos <- read.vcf.genotypes(filename)
#' 
#' # Display matrix dimensions
#' print(dim(genos))
#' 
#' # Display first 5 individuals and first 5 variants
#' print(genos[1:5, 1:5])
#'
#' @export
read.vcf.genotypes <- function(filename, regions) {
  if(missing(regions)) 
    regions <- character(0)

  fi <- path.expand(filename)
  if(!file.exists(fi))
    stop("File doesn't exist")

  readVCFgenotypes2(filename, regions)
}
