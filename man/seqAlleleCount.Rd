\name{seqAlleleCount}
\alias{seqAlleleCount}
\title{Get Allele Counts}
\description{
    Calculates the allele counts.
}
\usage{
seqAlleleCount(gdsfile, .progress=FALSE,
    parallel=getOption("seqarray.parallel", FALSE))
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{.progress}{if TRUE, show progress information}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (parallel
        processing) or other value; \code{parallel} is passed to the argument
        \code{cl} in \code{\link{seqParallel}}, see \code{\link{seqParallel}}
        for more details.}
}
\value{
    A list.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAlleleCount}}, \code{\link{seqNumAllele}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
f <- seqOpen(gds.fn)

head(seqAlleleCount(f, .progress=TRUE))


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
