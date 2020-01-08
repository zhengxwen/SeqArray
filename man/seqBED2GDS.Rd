\name{seqBED2GDS}
\alias{seqBED2GDS}
\alias{seqGDS2BED}
\title{Conversion between PLINK BED and SeqArray GDS}
\description{
    Conversion between PLINK BED format and SeqArray GDS format.
}
\usage{
seqBED2GDS(bed.fn, fam.fn, bim.fn, out.gdsfn, compress.geno="LZMA_RA",
    compress.annotation="LZMA_RA", optimize=TRUE, digest=TRUE, parallel=FALSE,
    verbose=TRUE)
seqGDS2BED(gdsfile, out.fn, multi.row=FALSE, verbose=TRUE)
}
\arguments{
    \item{bed.fn}{the file name of binary file, genotype information}
    \item{fam.fn}{the file name of first six columns of \code{".ped"}}
    \item{bim.fn}{the file name of extended MAP file:
        two extra columns = allele names}
    \item{gdsfile}{character (a GDS file name), or
        a \code{\link{SeqVarGDSClass}} object}
    \item{out.gdsfn}{the file name, output a file of SeqArray format}
    \item{out.fn}{the file name of PLINK binary format without extended names}
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
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (parallel
        processing), a numeric value indicating the number of cores, or a
        cluster object for parallel processing; \code{parallel} is passed to
        the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details}
    \item{multi.row}{if \code{TRUE}, a multiallelic site is converted to
        multiple rows in PLINK bim and bed files}
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

# convert bed to gds
seqBED2GDS(bed.fn, fam.fn, bim.fn, "tmp.gds")

seqSummary("tmp.gds")


# convert gds to bed
gdsfn <- seqExampleFileName("gds")
seqGDS2BED(gdsfn, "plink")


# remove the temporary file
unlink(c("tmp.gds", "plink.fam", "plink.bim", "plink.bed"), force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
