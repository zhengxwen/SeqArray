\name{seqExampleFileName}
\alias{seqExampleFileName}
\title{Example files}
\description{
    The example files of VCF and GDS format.
}
\usage{
seqExampleFileName(type=c("gds", "vcf", "KG_Phase1", "dosage"))
}
\arguments{
    \item{type}{"gds" (by default), "vcf", "KG_Phase1" or "dosage"}
}
\value{
    Return the path of a VCF file in the package if \code{type="vcf"}, or
the path of a GDS file if \code{type="gds"}; if \code{type="KG_Phase1"},
return the path of GDS file on Chromosome 22 of the 1000 Genomes Phase 1
project; if \code{type="dosage"}, return the path of GDS file with imputed
allele dosages.
}

\author{Xiuwen Zheng}

\examples{
seqExampleFileName("gds")

seqExampleFileName("vcf")

seqExampleFileName("KG_Phase1")

seqExampleFileName("dosage")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
