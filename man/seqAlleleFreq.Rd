\name{seqAlleleFreq}
\alias{seqAlleleFreq}
\title{Allele frequencies}
\description{
    Calculate the allele frequencies.
}
\usage{
seqAlleleFreq(gdsfile, ref.allele=0L)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{ref.allele}{\code{0} refers to the reference allele, \code{1} is the
        first alternative allele, ...; \code{NA}}
}
\value{
    A vector of allele frequencies.
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

summary(seqAlleleFreq(f, 0L))


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
