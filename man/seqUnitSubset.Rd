\name{seqUnitCreate}
\alias{seqUnitCreate}
\alias{seqUnitSubset}
\alias{seqUnitMerge}
\alias{seqUnitSetDiff}
\title{Subset and merge the units}
\description{
    Subset, merge and set-diff the variant unit(s).
}
\usage{
seqUnitCreate(idx, desp=NULL)
seqUnitSubset(units, i)
seqUnitMerge(...)
seqUnitSetDiff(ut1, ut2, drop=TRUE)
}
\arguments{
    \item{idx}{a list of numeric indexing vectors for specifying variants}
    \item{desp}{a \code{data.frame} for annotating the variant sets}
    \item{units}{a list of units of selected variants, with S3 class
        \code{SeqUnitListClass}}
    \item{...}{\code{SeqUnitListClass} objects to be merged}
    \item{i}{a numeric or logical vector for indices specifying elements}
    \item{ut1}{a \code{SeqUnitListClass} object}
    \item{ut2}{a \code{SeqUnitListClass} object, whose variant indices will
        be removed from \code{ut1}}
    \item{drop}{if \code{TRUE}, remove units that become empty after the
        set difference operation}
}
\value{
    The variant unit of \code{SeqUnitListClass}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqUnitSlidingWindows}}, \code{\link{seqUnitFilterCond}}
}

\examples{
# open the GDS file
gdsfile <- seqOpen(seqExampleFileName("gds"))

# variant units via sliding windows
units <- seqUnitSlidingWindows(gdsfile)

(u1 <- seqUnitSubset(units, 1:10))
(u2 <- seqUnitSubset(units, 30:39))

seqUnitMerge(u1, u2)

seqUnitCreate(list(1:10, 20:30), data.frame(gene=c("g1", "g2")))

# set difference: remove variants in u2 from u1
seqUnitSetDiff(u1, u2)

# close the GDS file
seqClose(gdsfile)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
