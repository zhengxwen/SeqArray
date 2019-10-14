\name{seqAddValue}
\alias{seqAddValue}
\title{Add values to a GDS File}
\description{
    Add or modify the values in a GDS file with hash code
}
\usage{
seqAddValue(gdsfile, varnm, val, desp="", replace=FALSE, compress="LZMA_RA",
    packed=FALSE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{varnm}{the variable name, e.g., "sample.id", "variant.id", "chromosome",
        "annotation/info/NEW_VARIABLE"}
    \item{val}{the R value can be integers, real numbers, characters,
        factor, logical or raw variable, \code{list} and \code{data.frame}}
    \item{desp}{variable description}
    \item{replace}{if \code{TRUE}, replace the existing variable silently if possible}
    \item{compress}{the compression method can be "" (no compression), see
        \code{\link{add.gdsn}}}
    \item{packed}{FALSE}
}
\value{
    Return none.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
library(SeqArray)

# the file of GDS
gds.fn <- seqExampleFileName("gds")
file.copy(gds.fn, "tmp.gds", overwrite=TRUE)

# display
(f <- seqOpen("tmp.gds", FALSE))

show(index.gdsn(f, "sample.id"))
seqAddValue(f, "sample.id", 1:90, replace=TRUE)
show(index.gdsn(f, "sample.id"))

show(index.gdsn(f, "chromosome"))
v <- seqGetData(f, "chromosome")
seqAddValue(f, "chromosome", paste0("chr", v), replace=TRUE)
show(index.gdsn(f, "chromosome"))
table(seqGetData(f, "chromosome"))

# sample annotation
seqAddValue(f, "sample.annotation", data.frame(ii=1:90, y=rep("A", 90)), replace=TRUE)
seqAddValue(f, "sample.annotation/float", (1:90)/90)

# close the GDS file
seqClose(f)


# remove the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
