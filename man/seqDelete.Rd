\name{seqDelete}
\alias{seqDelete}
\title{Delete GDS Variables}
\description{
    Deletes variables in the sequence GDS file.
}
\usage{
seqDelete(gdsfile, info.varname=character(), format.varname=character(),
    verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{info.varname}{the variables in the INFO field, i.e.,
        "annotation/info/VARIABLE_NAME"}
    \item{format.varname}{the variables in the FORMAT field, i.e.,
        "annotation/format/VARIABLE_NAME"}
    \item{verbose}{if \code{TRUE}, show information}
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
f

# close the GDS file
seqClose(f)

# clean up the fragments, reduce the file size
cleanup.gds("tmp.gds")


# remove the temporary file
unlink("tmp", force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
