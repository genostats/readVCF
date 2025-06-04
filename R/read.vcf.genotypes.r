
#' @export
read.vcf.genotypes <- function(filename, regions) {
  if(missing(regions)) 
    regions <- character(0)

  fi <- path.expand(filename)
  if(!file.exists(fi))
    stop("File doesn't exist")

  readVCFgenotypes2(filename, regions)
}
