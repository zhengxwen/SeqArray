\name{seqExport}
\alias{seqExport}
\title{Export to a GDS File}
\description{
    Exports to a GDS file with selected samples and variants, which are defined
by \code{seqSetFilter()}.
}
\usage{
seqExport(gdsfile, out.fn, info.var=NULL, fmt.var=NULL, samp.var=NULL,
    verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{out.fn}{the file name of output GDS file}
    \item{info.var}{characters, the variable name(s) in the INFO field
        for import; or \code{NULL} for all variables}
    \item{fmt.var}{characters, the variable name(s) in the FORMAT field
        for import; or \code{NULL} for all variables}
    \item{samp.var}{characters, the variable name(s) in the folder
        \code{"sample.annotation"}}
    \item{verbose}{if TRUE, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
# open the GDS file
(gds.fn <- seqExampleFileName("gds"))
(f <- seqOpen(gds.fn))

# get 'sample.id'
head(samp.id <- seqGetData(f, "sample.id"))

# get 'variant.id'
head(variant.id <- seqGetData(f, "variant.id"))

set.seed(100)
# set sample and variant filters
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8,10,12,14,16)])
seqSetFilter(f, variant.id=sample(variant.id, 100))


# export
seqExport(f, "tmp.gds")

(f1 <- seqOpen("tmp.gds")); seqClose(f1)


# close
seqClose(f)

# delete the temporary file
unlink("tmp.gds")
}

\keyword{gds}
\keyword{VCF}
\keyword{sequencing}
\keyword{genetics}
