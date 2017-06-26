\name{seqBlockApply}
\alias{seqBlockApply}
\title{Apply Functions Over Array Margins via Blocking}
\description{
    Returns a vector or list of values obtained by applying a function to
margins of genotypes and annotations via blocking.
}
\usage{
seqBlockApply(gdsfile, var.name, FUN, margin=c("by.variant"),
    as.is=c("none", "list", "unlist"), var.index=c("none", "relative", "absolute"),
    bsize=1024L, parallel=FALSE, .useraw=FALSE, .progress=FALSE, ...)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{var.name}{the variable name(s), see details}
    \item{FUN}{the function to be applied}
    \item{margin}{giving the dimension which the function will be applied over}
    \item{as.is}{returned value: a list, a vector or none;
        \code{as.is} can be a \code{\link{connection}} object, or a GDS node
        \code{\link{gdsn.class}} object; if "unlist" is used, produces a vector
        which contains all the atomic components, via
        \code{unlist(..., recursive=FALSE)}}
    \item{var.index}{if \code{"none"}, call \code{FUN(x, ...)} without variable
        index; if \code{"relative"} or \code{"absolute"}, add an argument to
        the user-defined function \code{FUN} like \code{FUN(index, x, ...)}
        where \code{index} is an index of variant starting from 1 if
        \code{margin = "by.variant"}: \code{"relative"} for indexing in the
        selection defined by \code{\link{seqSetFilter}}, \code{"absolute"}
        for indexing with respect to all data}
    \item{bsize}{block size}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details.}
    \item{.useraw}{\code{TRUE}, force to use RAW instead of INTEGER for
        genotypes and dosages; \code{FALSE}, use INTEGER; \code{NA}, use RAW
        instead of INTEGER if possible; for genotypes, 0xFF is missing value
        if RAW is used}
    \item{.progress}{if \code{TRUE}, show progress information}
    \item{...}{optional arguments to \code{FUN}}
}
\details{
    The variable name should be \code{"sample.id"}, \code{"variant.id"},
\code{"position"}, \code{"chromosome"}, \code{"allele"}, \code{"genotype"},
\code{"annotation/id"}, \code{"annotation/qual"}, \code{"annotation/filter"},
\code{"annotation/info/VARIABLE_NAME"}, or
\code{"annotation/format/VARIABLE_NAME"}.

\code{"$dosage"} is also allowed for the dosages of reference allele (integer:
0, 1, 2 and NA for diploid genotypes).

\code{"$num_allele"} returns an integer vector with the numbers of distinct
alleles.

    The algorithm is highly optimized by blocking the computations to exploit
the high-speed memory instead of disk.
}
\value{
    A vector, a list of values or none.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqApply}}, \code{\link{seqSetFilter}},
    \code{\link{seqGetData}}, \code{\link{seqParallel}},
    \code{\link{seqGetParallel}}
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
set.seed(100)
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8,10)],
    variant.id=sample(variant.id, 10))

# read
seqApply(f, "genotype", FUN=print, margin="by.variant")
seqApply(f, "genotype", FUN=print, margin="by.variant", .useraw=TRUE)

seqApply(f, "genotype", FUN=print, margin="by.sample")
seqApply(f, "genotype", FUN=print, margin="by.sample", .useraw=TRUE)

# read in block
seqGetData(f, "$dosage")
seqBlockApply(f, "$dosage", print, bsize=3)
seqBlockApply(f, "$dosage", function(x) x, as.is="list", bsize=3)
seqBlockApply(f, c(dos="$dosage", pos="position"), print, bsize=3)


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
