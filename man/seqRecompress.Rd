\name{seqRecompress}
\alias{seqRecompress}
\title{Recompress the GDS file}
\description{
    Recompress the SeqArray GDS file.
}
\usage{
seqRecompress(gds.fn, compress=c("ZIP", "LZ4", "LZMA", "Ultra", "UltraMax"),
    exclude=character(), optimize=TRUE, verbose=TRUE)
}
\arguments{
    \item{gds.fn}{the file name of SeqArray file}
    \item{compress}{the compression method}
    \item{exclude}{excluded GDS nodes}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    This function requires gdsfmt (>= 1.17.2).
}
\value{
    None.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}, \code{\link{seqStorageOption}}
}

\examples{
gds.fn <- seqExampleFileName("gds")

file.copy(gds.fn, "tmp.gds")

\dontrun{
seqRecompress(gds.fn, "UltraMax")
}

unlink("tmp.gds")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
