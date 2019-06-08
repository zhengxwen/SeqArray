\name{seqParallel}
\alias{seqParallel}
\alias{seqParApply}
\title{Apply Functions in Parallel}
\description{
    Applies a user-defined function in parallel.
}
\usage{
seqParallel(cl=seqGetParallel(), gdsfile, FUN,
    split=c("by.variant", "by.sample", "none"), .combine="unlist",
    .selection.flag=FALSE, .initialize=NULL, .finalize=NULL, .initparam=NULL,
    .balancing=FALSE, .bl_size=10000L, .bl_progress=FALSE, ...)
seqParApply(cl=seqGetParallel(), x, FUN, load.balancing=TRUE, ...)
}
\arguments{
    \item{cl}{\code{NULL} or \code{FALSE}: serial processing; \code{TRUE}:
        multicore processing (the maximum number of cores minor one);
        a numeric value: the number of cores to be used; a \code{cluster}
        object for parallel processing, created by the functions in the package
        \code{\link{parallel}}, like \code{\link{makeCluster}};
        a \code{BiocParallelParam} object from the BiocParallel package.
        See details}
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object, or \code{NULL}}
    \item{FUN}{the function to be applied, should be like
        \code{FUN(gdsfile, ...)} or \code{FUN(...)}}
    \item{split}{split the dataset by variant or sample according to multiple
        processes, or "none" for no split}
    \item{.combine}{define a fucntion for combining results from different
        processes; by default, \code{"unlist"} is used, to produce a vector
        which contains all the atomic components, via
        \code{unlist(..., recursive=FALSE)}; \code{"list"}, return a list
        of results created by child processes; \code{"none"}, no return; or a
        function, like "+".}
    \item{.selection.flag}{\code{TRUE} -- passes a logical vector of selection
        to the second argument of \code{FUN(gdsfile, selection, ...)}}
    \item{.initialize}{a user-defined function for initializing workers, should have
        two arguments (process_id, param)}
    \item{.finalize}{a user-defined function for finalizing workers, should have
        two arguments (process_id, param)}
    \item{.initparam}{parameters passed to \code{.initialize} and \code{.initialize}}
    \item{.balancing}{load balancing if \code{TRUE}}
    \item{.bl_size}{the increment for load balancing, 10000 for variants}
    \item{.bl_progress}{if \code{TRUE} and \code{.balancing=TRUE}, show progress
        information}
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
\code{SeqArray:::process_count} and \code{SeqArray:::process_index} to
tell the total number of cluster nodes and which cluster node being used.

    \code{seqParallel(, gdsfile=NULL, FUN=..., split="none")} could be used to
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
    \code{\link{seqApply}}, \code{\link{seqParallelSetup}},
    \code{\link{seqGetParallel}}
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
        paste(SeqArray:::process_index, SeqArray:::process_count, sep=" / ")
    }, split = "none")


################################################################


# close the GDS file
seqClose(gdsfile)


seqParallelSetup(FALSE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
