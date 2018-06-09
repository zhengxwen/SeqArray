\name{seqCheck}
\alias{seqCheck}
\title{Data Integrity Checking}
\description{
    Performs data integrity on a SeqArray GDS file.
}
\usage{
seqCheck(gdsfile, verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object, or a file name}
    \item{verbose}{if \code{TRUE}, display information}
}
\value{
    A list of the following components:
    \item{hash}{a \code{data.frame} for hash checking, including \code{algo}
        for digest algorithms and \code{ok} for the checking states}
    \item{dimension}{a \code{data.frame} for checking the dimension of each
        variable, including \code{ok} for the checking states and \code{info}
        for the error messages}
}

\author{Xiuwen Zheng}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

seqCheck(gds.fn)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
