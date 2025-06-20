#' Count variants in a VCF file
#'
#' @description
#' Efficiently counts the number of variant records in a VCF file, optionally limited to specific genomic regions.
#' This function supports both plain-text and bgzipped VCF files with an accompanying Tabix index (`.tbi`).
#'
#' @param filename A character string specifying the path to the VCF file (e.g., `"data/sample.vcf.gz"`).
#'                 The file may be compressed with bgzip and must have a corresponding `.tbi` index file if `regions` are used.
#' @param regions An optional character vector of genomic regions to restrict the count, formatted as `"chr:start-end"` 
#'                (e.g., `"chr2:10000-20000"`). If omitted, the entire file is scanned.
#'
#' @return A single integer representing the number of variants found in the file or specified regions.
#'
#' @details
#' When using the `regions` parameter, the regions must exactly match those indexed in the `.tbi` file. Use [getRegions()]
#' to inspect which regions are supported. Providing unsupported regions will result in an error such as:
#' `"Mismatch between regions and .tbi file ! Aborting"`.
#'
#' This function does not load variant content into memory, and is suitable for large-scale file scanning.
#'
#' @seealso [getRegions()], [openVCF()], [read.vcf.genotypes()]
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' countVariants(filename)
#'
#' @export
countVariants <- function(filename, regions) {
  if(missing(regions))
    regions <- character(0)

  fi <- path.expand(filename)
  if(!file.exists(fi))
    stop("File doesn't exist")

  countVariants_(fi, regions)
}

