
#' @export
getRegions <- function(x) {
  R <- getRegions_(x)
  chr <- getChroms(x)
  data.frame( chr = chr[R$tid + 1L], beg = R$beg, end = R$end )
}
