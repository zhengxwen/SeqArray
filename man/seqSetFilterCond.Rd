\name{seqSetFilterCond}
\alias{seqSetFilterCond}
\title{Set a Filter to Variant with Allele Count/Freq}
\description{
    Sets a filter to variant with specified allele count/frequency and
missing rate.
}
\usage{
seqSetFilterCond(gdsfile, maf=NaN, mac=1L, missing.rate=NaN,
    parallel=seqGetParallel(), .progress=FALSE, verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{maf}{minimum minor reference allele frequency}
    \item{mac}{minimum minor reference allele count}
    \item{missing.rate}{maximum missing genotype rate}
    \item{.progress}{if \code{TRUE}, show progress information}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details.}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    None.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilter}}, \code{\link{seqSetFilterChrom}},
    \code{\link{seqAlleleFreq}}, \code{\link{seqAlleleCount}},
    \code{\link{seqMissing}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
(f <- seqOpen(gds.fn))

seqSetFilterChrom(f, c(1, 6))
seqSetFilterCond(f, maf=0.05, .progress=TRUE)


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
