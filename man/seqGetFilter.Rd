\name{seqGetFilter}
\alias{seqGetFilter}
\title{Get the Filter of GDS File}
\description{
    Gets the filter of samples and variants.
}
\usage{
seqGetFilter(gdsfile, .useraw=FALSE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{.useraw}{returns logical vectors if \code{FALSE}, and returns
        raw vectors if \code{TRUE}}
}
\value{
    Return a list:
    \item{sample.sel}{a logical/raw vector indicating selected samples}
    \item{variant.sel}{a logical/raw vector indicating selected variants}
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilter}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

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


z <- seqGetFilter(f, .useraw=TRUE)
head(z$sample.sel)
head(z$variant.sel)


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
