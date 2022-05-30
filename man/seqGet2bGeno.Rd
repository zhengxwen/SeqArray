\name{seqGet2bGeno}
\alias{seqGet2bGeno}
\title{Get packed genotypes}
\description{
    Gets a RAW matrix of genotypes in a packed 2-bit format.
}
\usage{
seqGet2bGeno(gdsfile, samp_by_var=TRUE, ext_nbyte=0L, verbose=FALSE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{samp_by_var}{if \code{TRUE}, return a sample-by-variant matrix;
        otherwise, return a variant-by-sample matrix}
    \item{ext_nbyte}{additional \code{ext_nbyte} row(s) with missing genotypes}
    \item{verbose}{if \code{TRUE}, show progress information}
}
\value{
    Return a RAW matrix.
}
\details{
    If \code{samp_by_var=TRUE}, the function returns a sample-by-variant RAW
matrix (nrow = \code{ceiling(# of samples / 4)}); otherwise, it returns a
variant-by-sample RAW matrix (nrow = \code{ceiling(# of variants / 4)}). The
RAW matrix consists of a 2-bit array, with 0, 1 and 2 for dosage, and 3 for
missing genotype.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqGetData}}
}

\examples{
# open a GDS file
f <- seqOpen(seqExampleFileName("gds"))

str(seqGet2bGeno(f))

str(seqGet2bGeno(f, samp_by_var=FALSE))

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
