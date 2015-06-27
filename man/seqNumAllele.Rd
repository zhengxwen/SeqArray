\name{seqNumAllele}
\alias{seqNumAllele}
\title{Number of alleles}
\description{
    Return the numbers of alleles for each site.
}
\usage{
seqNumAllele(gdsfile)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
}
\value{
    The numbers of alleles for each site.
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
f <- seqOpen(gds.fn)

table(seqNumAllele(f))

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
