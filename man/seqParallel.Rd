\name{seqParallel}
\alias{seqParallel}
\title{Apply Functions in Parallel}
\description{
    Apply a user-defined function in parallel
}
\usage{
seqParallel(cl, gdsfile, FUN = function(gdsfile, ...) NULL,
    split=c("by.variant", "by.sample", "none"), .combine=NULL, ...)
}
\arguments{
    \item{cl}{NULL or a cluster object, created by the package \link[parallel]{parallel}
        or \href{http://CRAN.R-project.org/package=snow}{snow}}
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{FUN}{the function to be applied}
    \item{split}{split the dataset by variant or sample according to multiple
        processes, or "none" for no split}
    \item{.combine}{define a fucntion for combining results from different processes;
        by default, \code{'c'} is used; if \code{.combine==""}, return \code{invisible()}}
    \item{...}{optional arguments to \code{FUN}}
}
\details{
    If \code{cl = NULL} or \code{length(cl) == 0}, the function simply calls
\code{FUN(gdsfile, ...)}; otherwise, it splits jobs to different processes and calls
\code{FUN(gdsfile, ...)} on each process, the optional arguments are passed to different
processes.
}
\value{
    A vector or list of values.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilter}},
    \code{\link{seqGetData}}
    \code{\link{seqApply}}
}

\examples{
library(parallel)

# Use option cl.core to choose an appropriate cluster size or number of cores
cl <- makeCluster(getOption("cl.cores", 2))


# the file of GDS
gds.fn <- seqExampleFileName("gds")
# or gds.fn <- "C:/YourFolder/Your_GDS_File.gds"

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
