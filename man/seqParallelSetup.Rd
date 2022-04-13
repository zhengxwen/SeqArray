\name{seqParallelSetup}
\alias{seqParallelSetup}
\alias{seqGetParallel}
\alias{seqMulticoreSetup}
\title{Setup/Get a Parallel Environment}
\description{
    Setups a parallel environment in R for the current session.
}
\usage{
seqParallelSetup(cluster=TRUE, verbose=TRUE)
seqGetParallel()
seqMulticoreSetup(num, type=c("psock", "fork"), verbose=TRUE)
}
\arguments{
    \item{cluster}{\code{NULL} or \code{FALSE}: serial processing; \code{TRUE}:
        parallel processing with the maximum number of cores minor one;
        a numeric value: the number of cores to be used;
        a cluster object for parallel processing, created by the functions
        in the package \code{\link{parallel}}, like \code{\link{makeCluster}}.
        See details}
    \item{num}{the maximum number of cores used for the user-defined multicore
        setting; FALSE, NA or any value less than 2, to disable the multicore
        cluster}
    \item{type}{either PSOCK or Fork cluster setup for the multicore setting,
        the resulting parallel cluster will be used if 'parallel' is a number
        greater than one in associated functions}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    When \code{cl} is \code{TRUE} or a numeric value, forking techniques are
used to create a new child process as a copy of the current R process, see
\code{?parallel::mcfork}. However, forking is not available on Windows, so
multiple processes created by \code{\link{makeCluster}} are used instead.
    The R environment option \code{seqarray.parallel} will be set according
to the value of \code{cluster}.
    Using \code{seqParallelSetup(FALSE)} removes the registered cluster, as
does stopping the registered cluster.
}
\value{
    \code{seqParallelSetup()} has no return, and \code{seqGetParallel()} returns
\code{getOption("seqarray.parallel", FALSE)}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqParallel}}, \code{\link{seqApply}}
}

\examples{
library(parallel)

seqParallelSetup(2L)

# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
(f <- seqOpen(gds.fn))

# run in parallel
summary(seqMissing(f))

# close the GDS file
seqClose(f)

seqParallelSetup(FALSE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
