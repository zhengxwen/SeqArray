\name{seqDigest}
\alias{seqDigest}
\title{Hash function digests}
\description{
    Create hash function digests for all or a subset of data
}
\usage{
seqDigest(gdsfile, varname, algo=c("md5"), parallel=FALSE, verbose=FALSE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object or a GDS file name}
    \item{varname}{the variable name(s), see details}
    \item{algo}{the digest hash algorithm: "md5"}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details.}
    \item{verbose}{if \code{TRUE}, show progress information}
}
\details{
    The variable name should be \code{"sample.id"}, \code{"variant.id"},
\code{"position"}, \code{"chromosome"}, \code{"allele"}, \code{"annotation/id"},
\code{"annotation/qual"}, \code{"annotation/filter"},
\code{"annotation/info/VARIABLE_NAME"}, or
\code{"annotation/format/VARIABLE_NAME"}.

    Users can define a subset of data via \code{\link{seqSetFilter}} and
create a hash digest for the subset only.

    Note that the hashing algorithm is not parallelizable. When multiple
processes are used, the hash codes are calculated for each chuck of data, and
then call \code{digest()} to obtain the final hash code. Hence, the returned
hash codes are different when comparing \code{parallel=FALSE} and a multicore
implementation.
}
\value{
    A hash character.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilter}}, \code{\link{seqApply}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
f <- seqOpen(gds.fn)

seqDigest(f, "genotype")
seqDigest(f, "annotation/filter")
seqDigest(f, "annotation/format/DP")

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
