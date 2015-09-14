\name{seqParallel}
\alias{seqParallel}
\title{Apply Functions in Parallel}
\description{
    Applies a user-defined function in parallel.
}
\usage{
seqParallel(cl=getOption("seqarray.parallel", FALSE), gdsfile, FUN,
    split=c("by.variant", "by.sample", "none"), .combine="unlist",
    .selection.flag=FALSE, ...)
}
\arguments{
    \item{cl}{\code{NULL} or \code{FALSE}: serial processing; \code{TRUE}:
        parallel processing with the maximum number of cores minor one;
        a numeric value: the number of cores to be used;
        a cluster object for parallel processing, created by the functions
        in the package \code{\link{parallel}}, like \code{\link{makeCluster}}.
        See details}
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{FUN}{the function to be applied, should be like
        \code{FUN(gdsfile, ...)}}
    \item{split}{split the dataset by variant or sample according to multiple
        processes, or "none" for no split}
    \item{.combine}{define a fucntion for combining results from different
        processes; by default, \code{"unlist"} is used, to produce a vector
        which contains all the atomic components; \code{"list"}, return a list
        of results created by processes; \code{"none"}, no return; or a
        function, like "+".}
    \item{.selection.flag}{\code{TRUE} -- passes a logical vector of selection
        to the second argument of \code{FUN(gdsfile, selection, ...)}}
    \item{...}{optional arguments to \code{FUN}}
}
\details{
    When \code{cl} is \code{TRUE} or a numeric value, forking techniques are
used to create a new child process as a copy of the current R process, see
\code{?parallel::mcfork}. However, forking is not available on Windows, so
serial processing is used instead. In order to use multiple processes on
Windows, users have to create a cluster object via \code{\link{makeCluster}}.
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
library(parallel)

# Use option cl.core to choose an appropriate cluster size or number of cores
cl <- makeCluster(getOption("cl.cores", 2))


# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
(f <- seqOpen(gds.fn))

# the uniprocessor version
afreq1 <- seqParallel(NULL, f, FUN = function(gdsfile) {
        seqApply(gdsfile, "genotype", as.is="double",
            FUN=function(x) mean(x==0, na.rm=TRUE))
    }, split = "by.variant")

length(afreq1)
summary(afreq1)


# run in parallel
afreq2 <- seqParallel(cl, f, FUN = function(gdsfile) {
        seqApply(gdsfile, "genotype", as.is="double",
            FUN=function(x) mean(x==0, na.rm=TRUE))
    }, split = "by.variant")

length(afreq2)
summary(afreq2)


# check
all(afreq1 == afreq2)


################################################################
# check -- variant splits

seqParallel(cl, f, FUN = function(gdsfile) {
        v <- seqGetFilter(gdsfile)
        sum(v$variant.sel)
    }, split = "by.variant")
# [1] 674 674


################################################################

stopCluster(cl)

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
