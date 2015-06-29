\name{seqGetFilter}
\alias{seqGetFilter}
\title{Get the Filter of Samples and Variants}
\description{
    Get the filter of samples and variants.
}
\usage{
seqGetFilter(gdsfile)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
}
\value{
    Return a list:
    \item{sample.sel}{a logical vector for selected samples}
    \item{variant.sel}{a logical vector for selected variants}
}
\details{
    It is strongly suggested to call \code{seqOpen} instead of
\code{\link[gdsfmt]{openfn.gds}}, since \code{seqOpen} will initialize the
internal data for \code{seqGetData}, \code{seqApply}, etc.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqGetData}}, \code{\link{seqApply}}
}

\examples{
# the file of GDS
gds.fn <- seqExampleFileName("gds")
# or gds.fn <- "C:/YourFolder/Your_GDS_File.gds"

# display
(f <- seqOpen(gds.fn))

# get 'sample.id
(samp.id <- seqGetData(f, "sample.id"))
# "NA06984" "NA06985" "NA06986" ...

# get 'variant.id'
head(variant.id <- seqGetData(f, "variant.id"))


# set sample and variant filters
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8,10)])
set.seed(100)
seqSetFilter(f, variant.id=sample(variant.id, 10))

# get filter
z <- seqGetFilter(f)

# the number of selected samples
sum(z$sample.sel)
# the number of selected variants
sum(z$variant.sel)


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
