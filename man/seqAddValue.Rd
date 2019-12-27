\name{seqAddValue}
\alias{seqAddValue}
\title{Add values to a GDS File}
\description{
    Add or modify the values in a GDS file with hash code
}
\usage{
seqAddValue(gdsfile, varnm, val, desp=character(), replace=FALSE, compress="LZMA_RA",
    packed=TRUE, verbose=TRUE, verbose.attr=TRUE)
}
\arguments{
    \item{gdsfile}{character for file name, or a \code{\link{SeqVarGDSClass}} object}
    \item{varnm}{the variable name, e.g., "sample.id", "variant.id", "chromosome",
        "annotation/info/NEW_VARIABLE"}
    \item{val}{the R value can be integers, real numbers, characters,
        factor, logical or raw variable, \code{data.frame}}
    \item{desp}{variable description}
    \item{replace}{if \code{TRUE}, replace the existing variable silently if possible}
    \item{compress}{the compression method can be "" (no compression), see
        \code{\link{add.gdsn}}}
    \item{packed}{\code{TRUE}, pack data if there is any missing value}
    \item{verbose}{if \code{TRUE}, show information}
    \item{verbose.attr}{if \code{TRUE}, show attribute information in a GDS node}
}
\value{
    Return none.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
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

# annotation info
seqAddValue(f, "annotation/info/folder", NULL)  # add a new folder
seqAddValue(f, "annotation/info/folder/val", 1:1348, "random number")
seqAddValue(f, "annotation/info/folder/packed", c(rep(2L, 1000), rep(NA, 348)))
seqAddValue(f, "annotation/info/newff",
    data.frame(x=1:1348, y=rep("s", 1348), stringsAsFactors=FALSE),
    desp=c("integer numbers", "character"))

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
