\name{seqSlidingWindow}
\alias{seqSlidingWindow}
\title{Apply functions via a sliding window over variants}
\description{
    Returns a vector or list of values obtained by applying a function to a
sliding window over variants.
}
\usage{
seqSlidingWindow(gdsfile, var.name, win.size, shift=1, FUN,
    as.is = c("list", "integer", "double", "character", "none"),
    var.index = c("none", "relative", "absolute"), ...)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{var.name}{the variable name(s), see details}
    \item{win.size}{the size of sliding window}
    \item{shift}{the number of variants to shift the window at each step}
    \item{FUN}{the function to be applied}
    \item{as.is}{returned value: a list, an integer vector, etc}
    \item{var.index}{if \code{"none"}, call \code{FUN(x, ...)} without variable
        index; if \code{"relative"} or \code{"absolute"}, add an argument to
        the user-defined function \code{FUN} like \code{FUN(index, x, ...)}
        where \code{index} is an index of variant starting from 1:
        \code{"relative"} for indexing in the selection defined by
        \code{\link{seqSetFilter}}, \code{"absolute"} for indexing with
        respect to all data}
    \item{...}{optional arguments to \code{FUN}}
}
\details{
    The variable name should be \code{"sample.id"}, \code{"variant.id"},
\code{"position"}, \code{"chromosome"}, \code{"allele"}, \code{"annotation/id"},
\code{"annotation/qual"}, \code{"annotation/filter"},
\code{"annotation/info/VARIABLE_NAME"}, or \code{"annotation/format/VARIABLE_NAME"}.

    In the user-defined funciton \code{FUN(x, ...)} or
\code{FUN(index, x, ...)}, \code{x} is a list with \code{win.size} elements,
and each element includes values for the variable(s) \code{var.name};
\code{index} is the starting position of the sliding window.

    The algorithm is highly optimized by blocking the computations to exploit the
high-speed memory instead of disk.
}
\value{
    A vector or list of values.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilter}}, \code{\link{seqGetData}},
    \code{\link{seqApply}}
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
set.seed(100)
seqSetFilter(f, sample.id=samp.id[seq(2, 16, 2)],
    variant.id=sample(variant.id, 10))


# apply a function via a sliding window over variants
seqSlidingWindow(f, c(qual="annotation/id"), win.size=3,
    FUN = function(x) {
        # x is a list with 'win.size' elements
        print(x)
    }, as.is="none")

# apply a function via a sliding window over variants
seqSlidingWindow(f, c(qual="annotation/id"), win.size=3,
    FUN = function(x) {
        cat(unlist(x), sep="\t"); cat("\n")
    }, as.is="none")


# apply a function via a sliding window over variants
seqSlidingWindow(f, c(geno="genotype", phase="phase", qual="annotation/id"),
    FUN = function(index, x) {
        cat("Window ", index, ":\n", sep="")
        print(x)
    },
    win.size=3, as.is="none", var.index="relative")


# apply a function via a sliding window over variants
seqSlidingWindow(f, "genotype", win.size=4,
    FUN = function(index, x) {
        z <- unlist(lapply(x, function(z) mean(z, na.rm=TRUE)))
        cat("Window ", index, ", starting from Variant ", index,
            "\n    ", format(round(z,3), nsmall=3, width=8), "\n", sep="")
    },
    as.is="none", var.index="relative")


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
