#' @export
countVariants <- function(filename, regions) {
  if(missing(regions))
    regions <- character(0)

  fi <- path.expand(filename)
  if(!file.exists(fi))
    stop("File doesn't exist")

  countVariants_(fi, regions)
}

