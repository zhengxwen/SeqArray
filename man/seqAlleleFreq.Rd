\name{seqAlleleFreq}
\alias{seqAlleleFreq}
\alias{seqAlleleCount}
\title{Get Allele Frequencies or Counts}
\description{
    Calculates the allele frequencies or counts.
}
\usage{
seqAlleleFreq(gdsfile, ref.allele=0L, .progress=FALSE, parallel=seqGetParallel())
seqAlleleCount(gdsfile, ref.allele=0L, .progress=FALSE, parallel=seqGetParallel())
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{ref.allele}{\code{NULL}, a single numeric value, a numeric vector
        or a character vector; see Value}
    \item{.progress}{if \code{TRUE}, show progress information}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details.}
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
head(seqAlleleFreq(f, NULL, .progress=TRUE))

# return a numeric vector
summary(seqAlleleFreq(f, 0L, .progress=TRUE))

# return a numeric vector, AA is ancestral allele
AA <- seqGetData(f, "annotation/info/AA", .padNA=TRUE)
head(seqAlleleFreq(f, AA))

# allele counts
head(seqAlleleCount(f, NULL, .progress=TRUE))
head(seqAlleleCount(f, 0L, .progress=TRUE))
head(seqAlleleCount(f, AA, .progress=TRUE))

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
