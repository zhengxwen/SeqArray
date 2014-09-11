\name{seqExampleFileName}
\alias{seqExampleFileName}
\title{Example GDS file}
\description{
	The example GDS file for sequencing variants
}
\usage{
seqExampleFileName(type=c("gds", "vcf"))
}
\arguments{
	\item{type}{either "gds" or "vcf"}
}
\value{
	Return the file name of a VCF file shipped with the package if \code{type = "vcf"},
or the file name of a GDS file if \code{type = "gds"}.
}
\details{
	A GDS sequencing-variant file was created from a subset of VCF data of
the 1000 Genomes Project.
}

\author{Xiuwen Zheng}

\examples{
seqExampleFileName("gds")

seqExampleFileName("vcf")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
