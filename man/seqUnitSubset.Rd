\name{seqUnitCreate}
\alias{seqUnitCreate}
\alias{seqUnitSubset}
\alias{seqUnitMerge}
\title{Subset and merge the units}
\description{
    Subset and merge the variant unit(s).
}
\usage{
seqUnitCreate(idx, desp)
seqUnitSubset(units, i)
seqUnitMerge(ut1, ut2)
}
\arguments{
    \item{idx}{a list of numeric indexing vectors for specifying variants}
    \item{desp}{a \code{data.frame} for annotating the variant sets}
    \item{units}{a list of units of selected variants, with S3 class
        \code{SeqUnitListClass}}
    \item{ut1}{a list of units of selected variants, with S3 class
        \code{SeqUnitListClass}}
    \item{ut2}{a list of units of selected variants, with S3 class
        \code{SeqUnitListClass}}
    \item{i}{a numeric or logical vector for indices specifying elements}
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

# close the GDS file
seqClose(gdsfile)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
