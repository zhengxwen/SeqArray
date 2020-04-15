\name{seqAlleleFreq}
\alias{seqAlleleFreq}
\alias{seqAlleleCount}
\title{Get Allele Frequencies or Counts}
\description{
    Calculates the allele frequencies or counts for reference or minor alleles.
}
\usage{
seqAlleleFreq(gdsfile, ref.allele=0L, minor=FALSE, .progress=FALSE,
    parallel=seqGetParallel(), verbose=FALSE)
seqAlleleCount(gdsfile, ref.allele=0L, minor=FALSE, .progress=FALSE,
    parallel=seqGetParallel(), verbose=FALSE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{ref.allele}{\code{NULL}, a single numeric value, a numeric vector
        or a character vector; see Value}
    \item{minor}{if \code{TRUE}, return minor allele frequency/count}
    \item{.progress}{[deprecated] if TRUE, show progress information}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details.}
    \item{verbose}{if \code{TRUE}, show progress information}
}
\value{
    If \code{ref.allele=NULL}, the function returns a list of allele
frequencies/counts according to all allele per site. If \code{ref.allele} is a
single numeric value (like \code{0L}), it returns a numeric/integer vector for
the specified allele (\code{0L} for the reference allele, \code{1L} for the
first alternative allele, etc). If \code{ref.allele} is a numeric vector,
\code{ref.allele} specifies each allele per site. If \code{ref.allele} is a
character vector, \code{ref.allele} specifies the desired allele for each site
(e.g, ancestral allele for the derived allele frequency/count).
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqNumAllele}}, \code{\link{seqMissing}},
    \code{\link{seqParallel}}, \code{\link{seqGetParallel}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
f <- seqOpen(gds.fn)

# return a list
head(seqAlleleFreq(f, NULL, verbose=TRUE))

# return a numeric vector
summary(seqAlleleFreq(f, 0L, verbose=TRUE))

# return a numeric vector
summary(seqAlleleFreq(f, 0L, minor=TRUE, verbose=TRUE))

# return a numeric vector, AA is ancestral allele
AA <- seqGetData(f, "annotation/info/AA", .padNA=TRUE)
summary(seqAlleleFreq(f, AA))
summary(seqAlleleFreq(f, AA, minor=TRUE))

# allele counts
head(seqAlleleCount(f, NULL, verbose=TRUE))
head(seqAlleleCount(f, 0L, verbose=TRUE))
head(seqAlleleCount(f, 0L, minor=TRUE, verbose=TRUE))
head(seqAlleleCount(f, AA, verbose=TRUE))
head(seqAlleleCount(f, AA, minor=TRUE, verbose=TRUE))

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
