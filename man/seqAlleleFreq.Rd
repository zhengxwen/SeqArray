\name{seqAlleleFreq}
\alias{seqAlleleFreq}
\title{Allele frequencies}
\description{
    Calculate the allele frequencies.
}
\usage{
seqAlleleFreq(gdsfile, ref.allele=0L, parallel=FALSE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{ref.allele}{\code{NULL}, a single numeric value, a numeric vector
        or a character vector; see Value}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (parallel
        processing) or other value; \code{parallel} is passed to the argument
        \code{cl} in \code{\link{seqParallel}}, see \code{\link{seqParallel}}
        for more details.}
}
\value{
    If \code{ref.allele=NULL}, the function returns a list of allele
frequencies according to all allele per site. If \code{ref.allele} is a single
numeric value (like \code{0L}), it returns a numeric vector for the specified
alleles (\code{0L} for the reference allele, \code{1L} for the first alternative
allele, etc). If \code{ref.allele} is a numeric vector, \code{ref.allele}
specifies each allele per site. If \code{ref.allele} is a character vector,\
\code{ref.allele} specifies the desired allele for each site.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqNumAllele}}, \code{\link{seqMissing}},
    \code{\link{seqParallel}}
}

\examples{
# the file of GDS
gds.fn <- seqExampleFileName("gds")
# or gds.fn <- "C:/YourFolder/Your_GDS_File.gds"

# display
f <- seqOpen(gds.fn)

head(seqAlleleFreq(f, NULL))

summary(seqAlleleFreq(f, 0L))


AA <- toupper(seqGetData(f, "annotation/info/AA")$data)
summary(seqAlleleFreq(f, AA))


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
