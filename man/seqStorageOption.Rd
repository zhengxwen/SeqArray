\name{seqStorageOption}
\alias{seqStorageOption}
\title{Storage and Compression Options}
\description{
    Storage and compression options for GDS import and merging.
}
\usage{
seqStorageOption(compression=c("ZIP_RA.default", "ZIP_RA", "ZIP_RA.fast",
    "ZIP_RA.max", "LZ4_RA", "LZ4_RA.fast", "LZ4_RA.max", "none"),
    mode=NULL, float.mode="float32", geno.compress=NULL, info.compress=NULL,
    format.compress=NULL, index.compress=NULL, ...)
}
\arguments{
    \item{compression}{the default compression level, see \link{add.gdsn} for
        the description of compression methods; "ZIP_RA.default" is equivalent
        to "ZIP_RA"}
    \item{mode}{specify storage type for corresponding variable,
        e.g., `annotation/info/HM`="int16"}
    \item{float.mode}{specify the storage mode for read numbers, e.g.,
        "float32", "float64", "packedreal16"; the additional parameters can
        follow by colon, like "packedreal16:scale=0.0001"}
    \item{geno.compress}{NULL for the default value, or the compression method
        for genotypic data}
    \item{info.compress}{NULL for the default value, or the compression method
        for data sets stored in the INFO field (i.e., "annotation/info")}
    \item{format.compress}{NULL for the default value, or the compression method
        for data sets stored in the FORMAT field (i.e., "annotation/format")}
    \item{index.compress}{NULL for the default value, or the compression method
        for data index variables (e.g., "annotation/info/@HM")}
    \item{...}{other specified storage compression for corresponding variable,
        e.g., `annotation/info/HM`="ZIP_MAX"}
}
\value{
    Return a list with a class name "SeqGDSStorageClass".
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}, \code{\link{seqMerge}}
}

\examples{
# the file of VCF
(vcf.fn <- seqExampleFileName("vcf"))

# convert
seqVCF2GDS(vcf.fn, "tmp1.gds", storage.option=seqStorageOption())
(f1 <- seqOpen("tmp1.gds"))

# convert (maximize the compression ratio)
seqVCF2GDS(vcf.fn, "tmp2.gds", storage.option=seqStorageOption("ZIP_RA.max"))
(f2 <- seqOpen("tmp2.gds"))

# does not compress the genotypic data
seqVCF2GDS(vcf.fn, "tmp3.gds", storage.option=
    seqStorageOption("ZIP_RA", geno.compress=""))
(f3 <- seqOpen("tmp3.gds"))

# compress with LZ4
seqVCF2GDS(vcf.fn, "tmp4.gds", storage.option=seqStorageOption("LZ4_RA"))
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
