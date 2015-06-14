\name{seqCompress.Option}
\alias{seqCompress.Option}
\title{Compression Options for Importing VCF File(s)}
\description{
	Get compression options for importing VCF file(s)
}
\usage{
seqCompress.Option(default="ZIP_RA.MAX", ...)
}
\arguments{
	\item{default}{the default compression level, see \link{add.gdsn} for
	    more options}
	\item{...}{optional arguments}
}
\value{
	Return a list with a class name "seqCompress.class".
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
seqVCF2GDS(vcf.fn, "tmp1.gds", compress.option = seqCompress.Option())
(f1 <- seqOpen("tmp1.gds"))


# does not compress the genotypic data
seqVCF2GDS(vcf.fn, "tmp2.gds", compress.option = seqCompress.Option(genotype=""))
(f2 <- seqOpen("tmp2.gds"))


# close and remove
seqClose(f1)
seqClose(f2)
unlink(c("tmp1.gds", "tmp2.gds"))
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
