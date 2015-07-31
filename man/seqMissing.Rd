\name{seqMissing}
\alias{seqMissing}
\title{Missing genotype percentage}
\description{
    Calculate the missing rates per variant or per sample.
}
\usage{
seqMissing(gdsfile, per.variant=TRUE,
    parallel=getOption("seqarray.parallel", FALSE))
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{per.variant}{missing rate per variant if \code{TRUE}, or
        missing rate per sample if \code{FALSE}}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (parallel
        processing) or other value; \code{parallel} is passed to the argument
        \code{cl} in \code{\link{seqParallel}}, see \code{\link{seqParallel}}
        for more details.}
}
\value{
    A vector of missing rates.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAlleleFreq}}, \code{\link{seqNumAllele}},
    \code{\link{seqParallel}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
(f <- seqOpen(gds.fn))

summary(seqMissing(f, TRUE))

summary(seqMissing(f, FALSE))

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
