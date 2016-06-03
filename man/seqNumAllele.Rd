\name{seqNumAllele}
\alias{seqNumAllele}
\title{Number of alleles}
\description{
    Returns the numbers of alleles for each site.
}
\usage{
seqNumAllele(gdsfile)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
}
\value{
    The numbers of alleles for each site.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAlleleFreq}}, \code{\link{seqMissing}},
    \code{\link{seqParallel}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
f <- seqOpen(gds.fn)

table(seqNumAllele(f))

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
