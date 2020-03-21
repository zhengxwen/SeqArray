\name{seqSNP2GDS}
\alias{seqSNP2GDS}
\title{Convert SNPRelate Format to SeqArray Format}
\description{
    Converts a SNP GDS file to a SeqArray GDS file.
}
\usage{
seqSNP2GDS(gds.fn, out.fn, storage.option="LZMA_RA", major.ref=TRUE,
    ds.type=c("packedreal16", "float", "double"), optimize=TRUE, digest=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{gds.fn}{the file name of SNP format}
    \item{out.fn}{the file name, output a file of SeqArray format}
    \item{storage.option}{specify the storage and compression options,
        \code{"LZMA_RA"} to use LZMA compression algorithm with higher
        compression ratio compared to \code{"ZIP_RA"}}
    \item{major.ref}{if TRUE, use the major allele as a reference allele;
        otherwise, use A allele in SNP GDS file as a reference allele}
    \item{ds.type}{applicable when import dosages, the data type for storing
        dosages; see \code{\link{add.gdsn}}; \code{ds.type="packedreal16"} by
        default}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{digest}{a logical value (TRUE/FALSE) or a character ("md5", "sha1",
        "sha256", "sha384" or "sha512"); add hash codes to the GDS file if TRUE
        or a digest algorithm is specified}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the file name of SeqArray file with an absolute path. If the input
file is genotype dosage, the dosage matrix is stored in the node
\code{annotation/format/DS} with the estimated dosage of alternative alleles.
Any value less than 0 or greater than 2 will be replaced by NaN.
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
