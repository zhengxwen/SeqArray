\name{seqStorage.Option}
\alias{seqStorage.Option}
\title{Storage and Compression Options for Importing VCF File(s)}
\description{
    Storage and Compression Options for Importing VCF File(s).
}
\usage{
seqStorage.Option(compression="ZIP_RA.MAX", float.mode="float32")
}
\arguments{
    \item{compression}{the default compression level, see \link{add.gdsn} for
        more options}
    \item{float.mode}{specify the storage mode for read numbers, e.g.,
        "float32", "float64", "packedreal16"}
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
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# convert
seqVCF2GDS(vcf.fn, "tmp1.gds", storage.option = seqStorage.Option())
(f1 <- seqOpen("tmp1.gds"))


# does not compress the genotypic data
seqVCF2GDS(vcf.fn, "tmp2.gds", storage.option = seqStorage.Option())
(f2 <- seqOpen("tmp2.gds"))


# close and remove
seqClose(f1)
seqClose(f2)
unlink(c("tmp1.gds", "tmp2.gds"))
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
