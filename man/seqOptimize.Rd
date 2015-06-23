\name{seqOptimize}
\alias{seqOptimize}
\title{Transpose Data Array}
\description{
    Transpose data array or matrix for possibly higher-speed access
}
\usage{
seqOptimize(gdsfn, target=c("SeqVarTools"), format.var=TRUE, cleanup=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{gdsfn}{a \code{\link{SeqVarGDSClass}} object}
    \item{target}{"SeqVarTools" -- optimize GDS file for the SeqVarTools
        package}
    \item{format.var}{a character vector for selected variable names,
        or \code{TRUE} for all variables, according to "annotation/format"}
    \item{cleanup}{call \code{link{cleanup.gds}} if \code{TRUE}}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    None.
}
\details{
    It is designed for possibly higher-speed access. More details will be
provided in the future version.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqGetData}}, \code{\link{seqApply}}
}

\examples{
# the file name of VCF
(vcf.fn <- seqExampleFileName("vcf"))
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# convert
seqVCF2GDS(vcf.fn, "tmp.gds")

# prepare data for the SeqVarTools package
seqOptimize("tmp.gds", target="SeqVarTools")


# list the structure of GDS variables
(f <- seqOpen("tmp.gds"))
# close
seqClose(f)

unlink("tmp.gds")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
