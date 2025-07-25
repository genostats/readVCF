% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\docType{package}
\name{readVCF-package}
\alias{readVCF-package}
\alias{readVCF}
\title{Tools for Reading and Processing VCF Files}
\description{
The `readVCF` package provides efficient tools for reading and processing Variant Call Format (VCF) files, 
particularly useful for large-scale genotype data. It uses C++ backends for performance and integrates with 
`bigmemory` for scalable matrix storage.
}
\details{
Key features:
\itemize{
  \item Line-by-line VCF reading via [openVCF()]
  \item Extraction of sample and FORMAT data
  \item Genotype matrix construction and export to file-backed formats
  \item Compatible output (`.fam`, `.bim`, `.dosf`)
}
}
\author{
\strong{Maintainer}: Hervé Perdry \email{herve.perdry@universite-paris-saclay.fr}

Authors:
\itemize{
  \item Juliette Meyniel
  \item Maud Miclo
}

}
