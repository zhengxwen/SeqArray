\name{seqMissing}
\alias{seqMissing}
\title{Missing genotype percentage}
\description{
    Calculate the missing rates per variant or per sample.
}
\usage{
seqMissing(gdsfile, per.variant=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{per.variant}{missing rate per variant if \code{TRUE}, or
        missing rate per sample if \code{FALSE}}
}
\value{
    A vector of missing rates.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqOpen}}
}

\examples{
# the file of GDS
gds.fn <- seqExampleFileName("gds")
# or gds.fn <- "C:/YourFolder/Your_GDS_File.gds"

# display
(f <- seqOpen(gds.fn))

summary(seqMissing(f, TRUE))

summary(seqMissing(f, FALSE))

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
