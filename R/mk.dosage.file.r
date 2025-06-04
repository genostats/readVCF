#' 
#' @param vcf.file
#' @param dosage.file
#'
#' @examples # exemple provisoire
#' mk.dosage.file("inst/extdata/dosages.vcf.gz", "/tmp/essai" )
#' # attach it as a bigmemory object
#' dosages <- bigmemory::attach.big.matrix("/tmp/essai.dosf.desc")
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
