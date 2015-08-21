\name{seqVCF.SampID}
\alias{seqVCF.SampID}
\title{Get the Sample IDs}
\description{
    Returns the sample IDs of a VCF file.
}
\usage{
seqVCF.SampID(vcf.fn)
}
\arguments{
    \item{vcf.fn}{the file name of VCF}
}

\references{
    The variant call format and VCFtools.
    Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA,
    Handsaker RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R;
    1000 Genomes Project Analysis Group.
    Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF.Header}}, \code{\link{seqVCF2GDS}}
}

\examples{
# the VCF file
(vcf.fn <- seqExampleFileName("vcf"))

# get sample id
seqVCF.SampID(vcf.fn)
}

\keyword{gds}
\keyword{VCF}
\keyword{sequencing}
\keyword{genetics}
