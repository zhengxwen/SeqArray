\name{seqStorage.Option}
\alias{seqStorage.Option}
\title{Storage and Compression Options for Importing VCF File(s)}
\description{
    Storage and Compression Options for Importing VCF File(s).
}
\usage{
seqStorage.Option(compression="ZIP_RA", format.data="ZIP_RA:8M",
    float.mode="float32")
}
\arguments{
    \item{compression}{the default compression level, see \link{add.gdsn} for
        more options, e.g., "ZIP_RA", "ZIP_RA.max", "LZ4_RA", etc}
    \item{format.data}{the default compression level for data sets stored
        in the FORMAT field (i.e., "annotation/format")}
    \item{float.mode}{specify the storage mode for read numbers, e.g.,
        "float32", "float64", "packedreal16"; the additional parameters can
        follow by colon, like "packedreal16:scale=0.0001"}
}
\value{
    Return a list with a class name "SeqGDSStorageClass".
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
# the file of VCF
(vcf.fn <- seqExampleFileName("vcf"))

# convert
seqVCF2GDS(vcf.fn, "tmp1.gds", storage.option=seqStorage.Option())
(f1 <- seqOpen("tmp1.gds"))

# convert (maximize the compression ratio)
seqVCF2GDS(vcf.fn, "tmp2.gds", storage.option=seqStorage.Option("ZIP_RA.max",
    "ZIP_RA.max:8M"))
(f2 <- seqOpen("tmp2.gds"))

# does not compress the genotypic data
seqVCF2GDS(vcf.fn, "tmp3.gds", storage.option=
    seqStorage.Option(c("ZIP_RA", genotype="")))
(f3 <- seqOpen("tmp3.gds"))

# compress the FORMAT data with LZ4
seqVCF2GDS(vcf.fn, "tmp4.gds", storage.option=
    seqStorage.Option(format.data="LZ4"))
(f4 <- seqOpen("tmp4.gds"))


# close and remove the files
seqClose(f1)
seqClose(f2)
seqClose(f3)
seqClose(f4)

unlink(c("tmp1.gds", "tmp2.gds", "tmp3.gds", "tmp4.gds"))
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
