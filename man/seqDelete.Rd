\name{seqDelete}
\alias{seqDelete}
\title{Delete variables in a Sequencing GDS File}
\description{
	Delete variables in a sequencing GDS file
}
\usage{
seqDelete(gdsfile, info.varname=NULL, format.varname=NULL)
}
\arguments{
	\item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
	\item{info.varname}{the variables in the INFO field, i.e.,
		"annotation/info/VARIABLE_NAME"}
	\item{format.varname}{the variables in the FORMAT field, i.e.,
		"annotation/format/VARIABLE_NAME"}
}
\value{
	None.
}

\author{Xiuwen Zheng}
\seealso{
	\code{\link{seqOpen}}, \code{\link{seqClose}}
}

\examples{
# the file of VCF
vcf.fn <- seqExampleFileName("vcf")
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# convert
seqVCF2GDS(vcf.fn, "tmp.gds")

# display
(f <- seqOpen("tmp.gds", FALSE))

seqDelete(f, info.varname=c("HM2", "AA"), format.varname="DP")

# close the GDS file
seqClose(f)

# clean up the fragments, reduce the file size
cleanup.gds("tmp.gds")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
