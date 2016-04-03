\name{seqParallel}
\alias{seqParallel}
\alias{seqParApply}
\title{Apply Functions in Parallel}
\description{
    Applies a user-defined function in parallel.
}
\usage{
seqParallel(cl=getOption("seqarray.parallel", FALSE), gdsfile, FUN,
    split=c("by.variant", "by.sample", "none"), .combine="unlist",
    .selection.flag=FALSE, ...)
seqParApply(cl=getOption("seqarray.parallel", FALSE), x, FUN,
    load.balancing=TRUE, ...)
}
\arguments{
    \item{cl}{\code{NULL} or \code{FALSE}: serial processing; \code{TRUE}:
        parallel processing with the maximum number of cores minor one;
        a numeric value: the number of cores to be used;
        a cluster object for parallel processing, created by the functions
        in the package \code{\link{parallel}}, like \code{\link{makeCluster}}.
        See details}
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object, or \code{NULL}}
    \item{FUN}{the function to be applied, should be like
        \code{FUN(gdsfile, ...)} or \code{FUN(...)}}
    \item{split}{split the dataset by variant or sample according to multiple
        processes, or "none" for no split}
    \item{.combine}{define a fucntion for combining results from different
        processes; by default, \code{"unlist"} is used, to produce a vector
        which contains all the atomic components; \code{"list"}, return a list
        of results created by processes; \code{"none"}, no return; or a
        function, like "+".}
    \item{.selection.flag}{\code{TRUE} -- passes a logical vector of selection
        to the second argument of \code{FUN(gdsfile, selection, ...)}}
    \item{x}{a vector (atomic or list), passed to \code{FUN}}
    \item{load.balancing}{if \code{TRUE}, call \code{\link{clusterApplyLB}}
        instead of \code{\link{clusterApply}}}
    \item{...}{optional arguments to \code{FUN}}
}
\details{
    When \code{cl} is \code{TRUE} or a numeric value, forking techniques are
used to create a new child process as a copy of the current R process, see
\code{?parallel::mcfork}. However, forking is not available on Windows, and
\code{\link{makeCluster}} is called to make a cluster which will be
deallocated after calling \code{FUN}.

    It is strongly suggested to use \code{seqParallel} together with
\code{seqParallelSetup}. \code{seqParallelSetup} could work around the problem
of forking on Windows, without allocating clusters frequently.

    The user-defined function could use two predefined variables
\code{.seq_process_count} and \code{.seq_process_index} to tell the total
number of cluster nodes and which cluster node being used.

    \code{seqParallel(, gdsfile=NULL, FUN=..., split="none")} might be used to
setup multiple streams of pseudo-random numbers, and see
\code{\link{nextRNGStream}} or \code{\link{nextRNGSubStream}} in the package
\code{parallel}.
}
\value{
    A vector or list of values.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilter}}, \code{\link{seqGetData}},
    \code{\link{seqApply}}, \code{\link{seqParallelSetup}}
}

\examples{
library(parallel)

# choose an appropriate cluster size or number of cores
seqParallelSetup(2)


# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
(gdsfile <- seqOpen(gds.fn))

# the uniprocessor version
afreq1 <- seqParallel(, gdsfile, FUN = function(f) {
        seqApply(f, "genotype", as.is="double",
            FUN=function(x) mean(x==0, na.rm=TRUE))
    }, split = "by.variant")

length(afreq1)
summary(afreq1)


# run in parallel
afreq2 <- seqParallel(, gdsfile, FUN = function(f) {
        seqApply(f, "genotype", as.is="double",
            FUN=function(x) mean(x==0, na.rm=TRUE))
    }, split = "by.variant")

length(afreq2)
summary(afreq2)


# check
length(afreq1)  # 1348
all(afreq1 == afreq2)

################################################################
# check -- variant splits

seqParallel(, gdsfile, FUN = function(f) {
        v <- seqGetFilter(f)
        sum(v$variant.sel)
    }, split = "by.variant")
# [1] 674 674


################################################################

seqParallel(, NULL, FUN = function() {
        paste(.seq_process_index, .seq_process_count, sep="/")
    }, split = "none")


################################################################


# close the GDS file
seqClose(gdsfile)


seqParallelSetup(FALSE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
