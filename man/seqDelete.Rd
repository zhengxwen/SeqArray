\name{seqDelete}
\alias{seqDelete}
\title{Delete GDS Variables}
\description{
    Deletes variables in the SeqArray GDS file.
}
\usage{
seqDelete(gdsfile, info.var=character(), fmt.var=character(),
    samp.var=character(), verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{info.var}{the variables in the INFO field, i.e.,
        "annotation/info/VARIABLE_NAME"}
    \item{fmt.var}{the variables in the FORMAT field, i.e.,
        "annotation/format/VARIABLE_NAME"}
    \item{samp.var}{the variables in the sample annotation field, i.e.,
        "sample.annotation/VARIABLE_NAME"}
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
# the file of GDS
gds.fn <- seqExampleFileName("gds")
file.copy(gds.fn, "tmp.gds", overwrite=TRUE)

# display
(f <- seqOpen("tmp.gds", FALSE))

seqDelete(f, info.var=c("HM2", "AA"), fmt.var="DP")
f

# close the GDS file
seqClose(f)

# clean up the fragments, reduce the file size
cleanup.gds("tmp.gds")


# remove the temporary file
unlink("tmp.gds", force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
