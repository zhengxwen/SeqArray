\name{seqExampleFileName}
\alias{seqExampleFileName}
\title{Example files}
\description{
    The example files of VCF and GDS format.
}
\usage{
seqExampleFileName(type=c("gds", "vcf", "KG_Phase1"))
}
\arguments{
    \item{type}{either "gds" or "vcf"}
}
\value{
    Return the file name of a VCF file shipped with the package if
\code{type = "vcf"}, or the file name of a GDS file if \code{type = "gds"}.
}
\details{
    The SeqArray GDS file was created from a subset of VCF data of the 1000
Genomes Phase 1 Project.
}

\author{Xiuwen Zheng}

\examples{
seqExampleFileName("gds")

seqExampleFileName("vcf")

seqExampleFileName("KG_Phase1")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
