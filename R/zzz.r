#' @useDynLib readVCF, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @importFrom utils write.table
NULL

#' @title read a line in a VCF file
#' @name getLine
#' @description Read current line in a readVCF object 
#' 
#' @param x description
#' @param field 
#' @details blabla (GT et DS)
#' 
#' @returns a list with components values, CHR ...
#' 
#' @examples a <- openVCF(filename, "19:320000-700000")
#' sample.ids <- getSamples(a)
#' geno <- matrix(nrow = 0, ncol = length(sample.ids), dimnames = list(character(0), sample.ids))
#' while( VCFnext(a) ) {
#'   line <- getLine(a, "DS")
#'   geno <- rbind(geno, line$values)
#' }
#' 
#' @export getLine
NULL

#' @export getSamples
NULL

#' @rdname getLine
#' @name VCFnext
#' @export VCFnext
NULL

#' @export getFormats
NULL

#' @export getChroms
NULL



