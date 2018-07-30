\name{seqGDS2SNP}
\alias{seqGDS2SNP}
\title{Convert to a SNP GDS File}
\description{
    Converts a SeqArray GDS file to a SNP GDS file.
}
\usage{
seqGDS2SNP(gdsfile, out.gdsfn, dosage=FALSE, compress.geno="LZMA_RA",
    compress.annotation="LZMA_RA", optimize=TRUE, verbose=TRUE)
}
\arguments{
    \item{gdsfile}{character (GDS file name), or
        a \code{\link{SeqVarGDSClass}} object}
    \item{out.gdsfn}{the file name, output a file of VCF format}
    \item{dosage}{a logical value, or characters for the variable name of
        dosage in the SeqArray file; if \code{FALSE} exports genotypes,
        otherwise exports dosages}
    \item{compress.geno}{the compression method for "genotype"; optional
        values are defined in the function \code{add.gdsn}}
    \item{compress.annotation}{the compression method for the GDS variables,
        except "genotype"; optional values are defined in the function
        \code{add.gdsn}}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{verbose}{if \code{TRUE}, show information}
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
    \code{\link{seqSNP2GDS}}, \code{\link{seqVCF2GDS}},
    \code{\link{seqGDS2VCF}}
}

\examples{
# the GDS file
gds.fn <- seqExampleFileName("gds")

seqGDS2SNP(gds.fn, "tmp.gds")


# delete the temporary file
unlink("tmp.gds")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
