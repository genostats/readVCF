#' Get indexed genomic regions from a VCFReader object
#'
#' @description
#' Retrieves the list of genomic regions that are effectively accessible from a VCF file via its Tabix index (`.tbi`).
#' This is especially useful when querying specific regions with [openVCF()], as only regions listed here can be requested.
#'
#' @param x A `VCFReader` object, typically created using [openVCF()].
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{chr}{Chromosome name (as character).}
#'   \item{beg}{Start position of the region (0-based, inclusive).}
#'   \item{end}{End position of the region (0-based, exclusive).}
#' }
#'
#' @details
#' The returned regions correspond to the genomic chunks indexed in the `.tbi` file associated with the VCF.
#' When using `openVCF(filename, regions = ...)`, the requested regions **must exactly overlap** one or more of the returned regions;
#' otherwise, an error such as `"Mismatch between regions and .tbi file ! Aborting"` will occur.
#'
#' @seealso \link{openVCF}, \link{getChroms}, \link{VCFReader}
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' reader <- openVCF(filename, regions = c("19:320000-700000"))
#' getRegions(reader)
#'
#' @export
getRegions <- function(x) {
  R <- getRegions_(x)
  chr <- getChroms(x)
  data.frame( chr = chr[R$tid + 1L], beg = R$beg, end = R$end )
}
