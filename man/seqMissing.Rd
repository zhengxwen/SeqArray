\name{seqMissing}
\alias{seqMissing}
\title{Missing genotype percentage}
\description{
    Calculates the missing rates per variant or per sample.
}
\usage{
seqMissing(gdsfile, per.variant=TRUE, .progress=FALSE,
    parallel=seqGetParallel(), verbose=FALSE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{per.variant}{missing rate per variant if \code{TRUE},
        missing rate per sample if \code{FALSE}, or calculating missing rates
        for variants and samples if \code{NA}}
    \item{.progress}{[deprecated] if TRUE, show progress information}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details.}
    \item{verbose}{if \code{TRUE}, show progress information}
}
\value{
    A vector of missing rates, or a \code{list(variant, sample)} for both
variants and samples.
}
\details{
    If the gds node 'genotype/data' (integer genotypes) is not available,
the node 'annotation/format/DS' (numeric genotype dosages for alternative
alleles) will be used to calculate allele frequencies. At a site, it assumes
'annotation/format/DS' stores the dosage of the 1st alternative allele in the
1st column, 2nd alt. allele in the 2nd column if it is multi-allelic, and so on.
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

summary(m1 <- seqMissing(f, TRUE, verbose=TRUE))
summary(m2 <- seqMissing(f, FALSE, verbose=TRUE))

str(m <- seqMissing(f, NA, verbose=TRUE))
identical(m1, m$variant)  # should be TRUE
identical(m2, m$sample)   # should be TRUE

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
