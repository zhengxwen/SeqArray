\name{seqUnitSlidingWindows}
\alias{seqUnitSlidingWindows}
\title{Sliding units of selected variants}
\description{
    Generating units of selected variants via sliding windows.
}
\usage{
seqUnitSlidingWindows(gdsfile, win.size=5000L, win.shift=2500L, win.start=0L,
    verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{win.size}{window size in basepair}
    \item{win.shift}{the shift of sliding window in basepair}
    \item{win.start}{the start position in basepair}
    \item{verbose}{if \code{TRUE}, display information}
}
\value{
    A list of integer vectors with the class name "SeqUnitListClass".
}

\author{Xiuwen Zheng}

\examples{
# open the GDS file
gdsfile <- seqOpen(seqExampleFileName("gds"))

v <- seqUnitSlidingWindows(gdsfile)

# close the GDS file
seqClose(gdsfile)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
