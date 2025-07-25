% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\name{getChroms}
\alias{getChroms}
\title{List chromosome names in the VCF file}
\arguments{
\item{x}{A `VCFReader` object.}
}
\value{
A character vector of chromosome names.
}
\description{
Returns the chromosome identifiers as defined in the VCF header.
}
\examples{
filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
a <- openVCF(filename)
getChroms(a)

}
\seealso{
\link{VCFReader}
}
