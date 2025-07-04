#' Convert VCF to file-backed dosage matrix 
#'
#' @description
#' Converts genotype data from a VCF file into a file-backed matrix (dosages), and generates associated .bim and .fam files.
#' This enables interoperability with tools, and allows access via the bigmemory package.
#' Dosages can be extracted from the `"DS"` (dosage) or `"GP"` (genotype probabilities) FORMAT field.
#'
#' @param vcf.file A character string specifying the path to the input VCF file (gzipped and tabix-indexed recommended).
#' @param dosage.file A string prefix for output files (will generate dosage.dosf, .bim, .fam, .dosf.desc).
#' @param regions An optional character vector of genomic regions to extract (e.g., "19:320000-700000"). If omitted, the full file is read.
#' @param type A string specifying the storage format for the dosage matrix. One of `"float"` (default, 32-bit float, `.dosf`) or `"integer"` (16-bit unsigned integer, `.dos16`).
#' @param field The FORMAT field used to extract dosage information. One of `"DS"` (default) or `"GP"`.
#'
#' @details
#' This function performs the following:
#' \itemize{
#'   \item Extracts genotype dosage values (typically from the "DS" FORMAT field).
#'   \item Stores them in a .dosf file as a binary float matrix.
#'   \item Creates a .bim file containing SNP metadata (e.g., ID, CHR, POS, A1, A2).
#'   \item Creates a .fam file with sample information (family ID, individual ID, father ID, mother ID, sex, phenotype).
#'   \item Writes a .desc file to allow attaching the .dosf file via [bigmemory::attach.big.matrix()].
#' }
#' Files that already exist are not overwritten.
#'
#' @return This function is called for its side effects. It returns invisible(NULL).
#'
#' @examples
#' outprefix_ds <- tempfile("dosage_ds_")
#' outprefix_gp <- tempfile("dosage_gp_")
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#'
#' # From DS field
#' mk.dosage.file(filename, outprefix_ds, regions = "19:320000-330000", field = "DS")
#'
#' # From GP field
#' mk.dosage.file(filename, outprefix_gp, regions = "19:320000-330000", field = "GP")
#'
#' # Attach with bigmemory
#' library(bigmemory)
#' dos_ds <- attach.big.matrix(paste0(outprefix_ds, ".dosf.desc"))
#' dos_gp <- attach.big.matrix(paste0(outprefix_gp, ".dosf.desc"))
#'
#' # Compare first few values
#' print(dos_ds[1:min(5, nrow(dos_ds)), 1:min(5, ncol(dos_ds))])
#' print(dos_gp[1:min(5, nrow(dos_gp)), 1:min(5, ncol(dos_gp))])
#' 
#' @seealso \link[bigmemory]{attach.big.matrix}, \link{mk.descriptor.file}
#'
#' @export
mk.dosage.file <- function(vcf.file, dosage.file, regions, type = c("float", "integer"), field = "DS") {
  if(missing(regions)) regions <- character(0)
  if (is.null(type)) type <- "float"
  type <- match.arg(type)  # sécurise l'argument
  field <- match.arg(field, choices = c("DS", "GP"))
  
  L <- switch(type,
              float   = writeDosage(vcf.file, dosage.file, regions, field),
              integer = writeDosageInteger(vcf.file, dosage.file, regions, field))
  
  
  
  # let's do dosage.file.fam
  fam <- data.frame(famid = L$samples, id = L$samples, father = 0, mother = 0, sex = 0, phen = 0)
  famfile <- paste0(dosage.file, ".fam")
  if(file.exists(famfile)) {
    warning(famfile, " already exists, won't erase")
  } else {
    write.table(fam, famfile, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # choose extension and type according to 'type'
  ext <- if (type == "float") ".dosf" else ".dos16"
  dtype <- if (type == "float") "float" else "short"
  
  # build the corresponding .desc
  mk.descriptor.file(paste0(dosage.file, ext), length(L$samples), L$nbSNPs, dtype)
  
}
