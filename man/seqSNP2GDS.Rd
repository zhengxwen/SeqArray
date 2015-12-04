\name{seqSNP2GDS}
\alias{seqSNP2GDS}
\title{Convert SNPRelate Format to SeqArray Format}
\description{
    Converts a SNP GDS file to a SeqArray GDS file.
}
\usage{
seqSNP2GDS(gds.fn, out.gdsfn, compress.geno="ZIP_RA",
    compress.annotation="ZIP_RA", optimize=TRUE, digest=TRUE, verbose=TRUE)
}
\arguments{
    \item{gds.fn}{the file name of SNP format}
    \item{out.gdsfn}{the file name, output a file of SeqArray format}
    \item{compress.geno}{the compression method for "genotype"; optional
        values are defined in the function \code{add.gdsn}}
    \item{compress.annotation}{the compression method for the GDS variables,
        except "genotype"; optional values are defined in the function
        \code{add.gdsn}}
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
