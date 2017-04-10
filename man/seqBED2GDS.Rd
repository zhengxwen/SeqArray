\name{seqBED2GDS}
\alias{seqBED2GDS}
\title{Convert PLINK BED Format to SeqArray Format}
\description{
    Converts a PLINK BED file to a SeqArray GDS file.
}
\usage{
seqBED2GDS(bed.fn, fam.fn, bim.fn, out.gdsfn,
    compress.geno="LZMA_RA", compress.annotation="LZMA_RA",
    optimize=TRUE, digest=TRUE, verbose=TRUE)
}
\arguments{
    \item{bed.fn}{the file name of binary file, genotype information}
    \item{fam.fn}{the file name of first six columns of \code{".ped"}}
    \item{bim.fn}{the file name of extended MAP file:
        two extra columns = allele names}
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
    \code{\link{seqSNP2GDS}}, \code{\link{seqVCF2GDS}}
}

\examples{
library(SNPRelate)

# PLINK BED files
bed.fn <- system.file("extdata", "plinkhapmap.bed.gz", package="SNPRelate")
fam.fn <- system.file("extdata", "plinkhapmap.fam.gz", package="SNPRelate")
bim.fn <- system.file("extdata", "plinkhapmap.bim.gz", package="SNPRelate")

# convert
seqBED2GDS(bed.fn, fam.fn, bim.fn, "tmp.gds")

seqSummary("tmp.gds")

# remove the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
