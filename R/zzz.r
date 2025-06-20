#' @docType package
#' @name readVCF-package
#' @aliases readVCF
#' @title Tools for Reading and Processing VCF Files
#' 
#' @description
#' The `readVCF` package provides efficient tools for reading and processing Variant Call Format (VCF) files, 
#' particularly useful for large-scale genotype data. It uses C++ backends for performance and integrates with 
#' `bigmemory` for scalable matrix storage.
#' 
#' @details
#' Key features:
#' \itemize{
#'   \item Line-by-line VCF reading via [openVCF()]
#'   \item Extraction of sample and FORMAT data
#'   \item Genotype matrix construction and export to file-backed formats
#'   \item Compatible output (`.fam`, `.bim`, `.dosf`)
#' }
#'
#' @useDynLib readVCF, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils write.table
"_PACKAGE"


#' Read the current variant line from a VCFReader
#'
#' @name getLine
#'
#' @description
#' Extracts the current line from a `VCFReader` object, returning the requested FORMAT field (e.g., `"GT"`, `"DS"`)
#' along with metadata such as chromosome and position.
#'
#' @param x A `VCFReader` object.
#' @param field A character string indicating the FORMAT field to extract (e.g., `"GT"`, `"DS"`).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{values}{A named vector of values (type depends on FORMAT field, e.g., integer for `"GT"`, numeric for `"DS"`). Empty if no samples or missing data.}
#'   \item{CHR}{Internal chromosome identifier (integer index as per tabix).}
#'   \item{POS}{Genomic position (1-based).}
#'   \item{ID, REF, ALT, QUAL, FILTER, INFO}{Character strings with standard VCF metadata fields. Empty strings if missing.}
#' }
#'
#' @details
#' Note that `CHR` is an internal numeric identifier from the tabix index.
#' To obtain chromosome names, use [getChroms()] with the index value:  
#' `getChroms(x)[ line$CHR + 1L ]`
#'
#' @seealso [VCFnext()], [getChroms()], [getSamples()], [getFormats()]
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' a <- openVCF(filename, "19:320000-700000")
#' VCFnext(a)
#' line <- getLine(a, "GT")
#' getChroms(a)[ line$CHR + 1L ]
#' line$values
#'
#' @export getLine
NULL


#' Get sample IDs from a VCFReader object
#'
#' @name getSamples
#'
#' @description
#' Extracts the list of sample identifiers from the header of a VCF file opened with [openVCF()].
#'
#' @param x A `VCFReader` object.
#'
#' @return A character vector of sample IDs.
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' a <- openVCF(filename)
#' getSamples(a)
#'
#' @export getSamples
NULL



#' Move to the next line in a VCFReader object
#'
#' @name VCFnext
#'
#' @description
#' Advances the internal pointer to the next line in the VCF file.
#' Should be used in a loop to iterate over lines.
#'
#' @param x A `VCFReader` object.
#'
#' @return `TRUE` if a new line was successfully read, `FALSE` if end of file is reached.
#'
#' @seealso [getLine()]
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' a <- openVCF(filename)
#' VCFnext(a)
#'
#' @export VCFnext
NULL


#' List available FORMAT fields in a VCFReader object
#'
#' @name getFormats
#'
#' @description
#' Returns the names of FORMAT fields available in the VCF header (e.g., "GT", "DS", "GP").
#'
#' @param x A `VCFReader` object.
#'
#' @return A character vector of FORMAT field names.
#'
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' a <- openVCF(filename)
#' getFormats(a)
#'
#' @export getFormats
NULL


#' List chromosome names in the VCF file
#'
#' @name getChroms
#'
#' @description
#' Returns the chromosome identifiers as defined in the VCF header.
#'
#' @param x A `VCFReader` object.
#'
#' @return A character vector of chromosome names.
#' 
#' @examples
#' filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
#' a <- openVCF(filename)
#' getChroms(a)
#'
#' @export getChroms
NULL



