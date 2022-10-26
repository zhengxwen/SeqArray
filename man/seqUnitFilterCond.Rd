\name{seqUnitFilterCond}
\alias{seqUnitFilterCond}
\title{Filter unit variants}
\description{
    Filters out the unit variants according to MAF, MAC and missing rates.
}
\usage{
seqUnitFilterCond(gdsfile, units, maf=NaN, mac=1L, missing.rate=NaN,
    minsize=1L, parallel=seqGetParallel(), verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{units}{a list of units of selected variants, with S3 class
        \code{SeqUnitListClass}}
    \item{maf}{minimum minor reference allele frequency, or a range of MAF
        \code{maf[1] <= ... < maf[2]}}
    \item{mac}{minimum minor reference allele count, or a range of MAC
        \code{mac[1] <= ... < mac[2]}}
    \item{missing.rate}{maximum missing genotype rate}
    \item{minsize}{the minimum of unit size}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details.}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    A S3 object with the class name "SeqUnitListClass" and two components
(\code{desp} and \code{index}): the first is a data.frame with columns "chr",
"start" and "end", and the second is list of integer vectors (the variant
indices).
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqUnitApply}}, \code{\link{seqUnitCreate}},
    \code{\link{seqUnitSubset}}, \code{\link{seqUnitMerge}}
}

\examples{
# open the GDS file
gdsfile <- seqOpen(seqExampleFileName("gds"))

unit1 <- seqUnitSlidingWindows(gdsfile)
unit1 # "desp"  "index"

# only rare variants
newunit <- seqUnitFilterCond(gdsfile, unit1, maf=c(0, 0.01))
newunit


# excluded variants
exvar <- setdiff(unique(unlist(unit1$index)), unique(unlist(newunit$index)))

seqSetFilter(gdsfile, variant.sel=exvar)
maf <- seqAlleleFreq(gdsfile, minor=TRUE)
table(maf > 0)
summary(maf[maf > 0])  # > 0.01


# close the GDS file
seqClose(gdsfile)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
