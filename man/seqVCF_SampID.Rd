\name{seqVCF_SampID}
\alias{seqVCF_SampID}
\title{Get the Sample IDs}
\description{
    Returns the sample IDs of a VCF file.
}
\usage{
seqVCF_SampID(vcf.fn)
}
\arguments{
    \item{vcf.fn}{the file name, output a file of VCF format; or a
        \code{\link{connection}} object}
}

\references{
    Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo,
M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., et al. (2011).
The variant call format and VCFtools. Bioinformatics 27, 2156-2158.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF_Header}}, \code{\link{seqVCF2GDS}}
}

\examples{
# the VCF file
(vcf.fn <- seqExampleFileName("vcf"))

# get sample id
seqVCF_SampID(vcf.fn)
}

\keyword{gds}
\keyword{VCF}
\keyword{sequencing}
\keyword{genetics}
