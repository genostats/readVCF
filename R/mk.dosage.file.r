#' Convert VCF to file-backed dosage matrix 
#'
#' @description
#' Converts genotype data from a VCF file into a file-backed matrix (dosages), and generates associated `.bim` and `.fam` files.
#' This enables interoperability with tools, and allows access via the `bigmemory` package.
#'
#' @param vcf.file A character string specifying the path to the input VCF file (gzipped and tabix-indexed recommended).
#' @param dosage.file A string prefix for output files (will generate `dosage.dosf`, `.bim`, `.fam`, `.dosf.desc`).
#' @param regions An optional character vector of genomic regions to extract (e.g., `"19:320000-700000"`). If omitted, the full file is read.
#'
#' @details
#' This function performs the following:
#' \itemize{
#'   \item Extracts genotype dosage values (typically from the `"DS"` FORMAT field).
#'   \item Stores them in a `.dosf` file as a binary float matrix.
#'   \item Creates a `.bim` file containing SNP metadata (e.g., ID, CHR, POS, A1, A2).
#'   \item Creates a `.fam` file with sample information (family ID, individual ID, father ID, mother ID, sex, phenotype).
#'   \item Writes a `.desc` file to allow attaching the `.dosf` file via [bigmemory::attach.big.matrix()].
#' }
#' Files that already exist are not overwritten.
#'
#' @return This function is called for its side effects. It returns `invisible(NULL)`.
#'
#' @examples
#' mk.dosage.file("inst/extdata/example.vcf.gz", "/tmp/essai")
#' # Attach it as a bigmemory matrix:
#' dosages <- bigmemory::attach.big.matrix("/tmp/essai.dosf.desc")
#'
#' @seealso [bigmemory::attach.big.matrix()], [mk.descriptor.file()]
#'
#' @export
mk.dosage.file <- function(vcf.file, dosage.file, regions) {
  if(missing(regions)) regions <- character(0)

  # this writes dosage.file.dosf / dosage.file.bim
  L <- writeDosage(vcf.file, dosage.file, regions)

  # let's do dosage.file.fam
  fam <- data.frame(famid = L$samples, id = L$samples, father = 0, mother = 0, sex = 0, phen = 0)
  famfile <- paste0(dosage.file, ".fam")
  if(file.exists(famfile)) {
    warning(famfile, " already exists, won't erase")
  } else {
    write.table(fam, famfile, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }

  # and descriptor...
  mk.descriptor.file(paste0(dosage.file, ".dosf"), length(L$samples), L$nbSNPs, "float")
}
