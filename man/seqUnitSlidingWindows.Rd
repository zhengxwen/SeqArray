\name{seqUnitSlidingWindows}
\alias{seqUnitSlidingWindows}
\title{Sliding units of selected variants}
\description{
    Generating units of selected variants via sliding windows.
}
\usage{
seqUnitSlidingWindows(gdsfile, win.size=5000L, win.shift=2500L, win.start=0L,
    dup.rm=TRUE, verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{win.size}{window size in basepair}
    \item{win.shift}{the shift of sliding window in basepair}
    \item{win.start}{the start position in basepair}
    \item{dup.rm}{if \code{TRUE}, remove duplicate and zero-length windows}
    \item{verbose}{if \code{TRUE}, display information}
}
\value{
    A S3 object with the class name "SeqUnitListClass" and two components
(\code{desp} and \code{index}): the first is a data.frame with columns "chr",
"start" and "end", and the second is list of integer vectors (the variant
indices).
}

\author{Xiuwen Zheng}

\examples{
# open the GDS file
gdsfile <- seqOpen(seqExampleFileName("gds"))

v <- seqUnitSlidingWindows(gdsfile)
names(v) # "desp"  "index"

head(v$desp)

# close the GDS file
seqClose(gdsfile)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
