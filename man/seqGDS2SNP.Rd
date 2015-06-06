\name{seqGDS2SNP}
\alias{seqGDS2SNP}
\title{Convert to a SNP GDS file}
\description{
	Convert a sequencing GDS file to a SNP GDS file
}
\usage{
seqGDS2SNP(gdsfile, out.gdsfn, compress.annotation="ZIP_RA.max",
    compress.geno="ZIP_RA.max", verbose=TRUE)
}
\arguments{
	\item{gdsfile}{character (GDS file name), or
	    a \code{\link{SeqVarGDSClass}} object}
	\item{out.gdsfn}{the file name, output a file of VCF format}
    \item{compress.annotation}{the compression method for the GDS variables,
        except "genotype"; optional values are defined in the function
        \code{add.gdsn}}
    \item{compress.geno}{the compression method for "genotype"; optional
        values are defined in the function \code{add.gdsn}}
	\item{verbose}{if TRUE, show information}
}
\value{
	Return the file name of VCF file with an absolute path.
}
\details{
	\code{\link{seqSetFilter}} can be used to define a subset of data for
the conversion.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}, \code{\link{seqGDS2VCF}}
}

\examples{
# the file of GDS
gds.fn <- seqExampleFileName("gds")

seqGDS2SNP(gds.fn, "test.gds")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
