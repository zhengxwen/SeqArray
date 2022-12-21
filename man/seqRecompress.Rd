\name{seqRecompress}
\alias{seqRecompress}
\title{Recompress the GDS file}
\description{
    Recompress the SeqArray GDS file.
}
\usage{
seqRecompress(gds.fn, compress=c("ZIP", "LZ4", "LZMA", "Ultra", "UltraMax", "none"),
    exclude=character(), optimize=TRUE, verbose=TRUE)
}
\arguments{
    \item{gds.fn}{the file name of SeqArray file}
    \item{compress}{the compression method, \code{compress="ZIP"} by default}
    \item{exclude}{a list of GDS nodes to be excluded, see details}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    This function requires gdsfmt (>= v1.17.2). \code{\link{seqVCF2GDS}} usually
takes lots of memory when the compression method \code{"LZMA_RA.max"},
\code{"Ultra"} or \code{"UltraMax"} is specified. So users could call
\code{seqVCF2GDS(, storage.option="ZIP_RA")} first, and then recompress the
GDS file with a higher compression option, e.g., \code{"UltraMax"}.
\code{seqRecompress()} takes much less memory than \code{seqVCF2GDS()},
since it recompresses data in a GDS node each time.

    \code{"UltraMax"} might be not better than \code{"Ultra"}, and its behavior
is similar to \code{xz -9 --extreme}: use a slower variant of the selected
compression preset level (-9) to hopefully get a little bit better compression
ratio, but with bad luck this can also make it worse.

    \code{ls.gdsn(gdsfile, include.hidden=TRUE, recursive=TRUE)} returns a list
of GDS nodes to be re-compressed, and users can specify the excluded nodes in
the argument \code{exclude}.
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

seqRecompress("tmp.gds", "LZMA")

unlink("tmp.gds")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
