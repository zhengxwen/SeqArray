\name{seqMissing}
\alias{seqMissing}
\title{Missing genotype percentage}
\description{
    Calculates the missing rates per variant or per sample.
}
\usage{
seqMissing(gdsfile, per.variant=TRUE, .progress=FALSE, parallel=seqGetParallel())
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{per.variant}{missing rate per variant if \code{TRUE}, or
        missing rate per sample if \code{FALSE}}
    \item{.progress}{if TRUE, show progress information}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details.}
}
\value{
    A vector of missing rates.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAlleleFreq}}, \code{\link{seqNumAllele}},
    \code{\link{seqParallel}}, \code{\link{seqGetParallel}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
(f <- seqOpen(gds.fn))

summary(seqMissing(f, TRUE, .progress=TRUE))

summary(seqMissing(f, FALSE, .progress=TRUE))

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
