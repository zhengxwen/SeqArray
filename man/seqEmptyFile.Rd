\name{seqEmptyFile}
\alias{seqEmptyFile}
\title{Empty GDS file}
\description{
    Create a new empty GDS file.
}
\usage{
seqEmptyFile(outfn, sample.id=character(), verbose=TRUE)
}
\arguments{
    \item{outfn}{the output file name for a GDS file}
    \item{sample.id}{a list of sample IDs}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    None.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
seqEmptyFile("tmp.gds")

(f <- seqOpen("tmp.gds"))
seqClose(f)

# remove the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
