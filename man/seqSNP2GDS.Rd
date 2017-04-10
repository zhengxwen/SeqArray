\name{seqSNP2GDS}
\alias{seqSNP2GDS}
\title{Convert SNPRelate Format to SeqArray Format}
\description{
    Converts a SNP GDS file to a SeqArray GDS file.
}
\usage{
seqSNP2GDS(gds.fn, out.fn, storage.option="LZMA_RA", optimize=TRUE,
    digest=TRUE, verbose=TRUE)
}
\arguments{
    \item{gds.fn}{the file name of SNP format}
    \item{out.fn}{the file name, output a file of SeqArray format}
    \item{storage.option}{specify the storage and compression options,
        \code{"LZMA_RA"} to use LZMA compression algorithm with higher
        compression ratio compared to \code{"ZIP_RA"}}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{digest}{a logical value (TRUE/FALSE) or a character ("md5", "sha1",
        "sha256", "sha384" or "sha512"); add hash codes to the GDS file if TRUE
        or a digest algorithm is specified}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the file name of SeqArray file with an absolute path.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqGDS2SNP}}, \code{\link{seqVCF2GDS}},
    \code{\link{seqGDS2VCF}}, \code{\link{seqBED2GDS}}
}

\examples{
library(SNPRelate)

# the GDS file
gds.fn <- snpgdsExampleFileName()

seqSNP2GDS(gds.fn, "tmp.gds")

seqSummary("tmp.gds")


# remove the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
