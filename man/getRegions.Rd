% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRegions.r
\name{getRegions}
\alias{getRegions}
\title{Get indexed genomic regions from a VCFReader object}
\usage{
getRegions(x)
}
\arguments{
\item{x}{A `VCFReader` object, typically created using [openVCF()].}
}
\value{
A data frame with three columns:
\describe{
  \item{chr}{Chromosome name (as character).}
  \item{beg}{Start position of the region (0-based, inclusive).}
  \item{end}{End position of the region (0-based, exclusive).}
}
}
\description{
Retrieves the list of genomic regions that are effectively accessible from a VCF file via its Tabix index (`.tbi`).
This is especially useful when querying specific regions with [openVCF()], as only regions listed here can be requested.
}
\details{
The returned regions correspond to the genomic chunks indexed in the `.tbi` file associated with the VCF.
When using `openVCF(filename, regions = ...)`, the requested regions **must exactly overlap** one or more of the returned regions;
otherwise, an error such as `"Mismatch between regions and .tbi file ! Aborting"` will occur.
}
\examples{
filename <- system.file("extdata", "example.vcf.gz", package = "readVCF")
reader <- openVCF(filename, regions = c("19:320000-700000"))
getRegions(reader)

}
\seealso{
\link{openVCF}, \link{getChroms}, \link{VCFReader}
}
