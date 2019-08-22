\name{seqAddValue}
\alias{seqAddValue}
\title{Add values to a GDS File}
\description{
    Add or modify the values in a GDS file with hash code
}
\usage{
seqAddValue(gdsfile, varnm, val, replace=FALSE, compress="LZMA_RA")
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{varnm}{the variable name, e.g., "sample.id", "variant.id", "chromosome",
        "annotation/info/new.variable"}
    \item{val}{the R value can be integers, real numbers, characters,
        factor, logical or raw variable, \code{list} and \code{data.frame}}
    \item{replace}{if \code{TRUE}, replace the existing variable silently if possible}
    \item{compress}{the compression method can be "" (no compression), see
        \code{\link{add.gdsn}}}
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


# close the GDS file
seqClose(f)


# remove the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
