\name{seqSetFilterChrom}
\alias{seqSetFilterChrom}
\title{Chromosome Selection}
\description{
    Selects the variants according to the specified chromosome(s).
}
\usage{
seqSetFilterChrom(gdsfile, include=NULL, is.num=NA, from.bp=NaN, to.bp=NaN,
    verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{include}{NULL, or character for specified chromosome(s)}
    \item{is.num}{a logical variable: \code{TRUE}, chromosome code is numeric;
        \code{FALSE}, chromosome is not numeric}
    \item{from.bp}{numeric, the lower bound of position}
    \item{to.bp}{numeric, the upper bound of position}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    None.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilter}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

f <- seqOpen(gds.fn)
seqSummary(f)


seqSetFilterChrom(f, is.num=TRUE)
seqSummary(f, "genotype", check="none")

seqSetFilterChrom(f, is.num=FALSE)
seqSummary(f, "genotype", check="none")

seqSetFilterChrom(f, 1:4)
seqSummary(f, "genotype", check="none")
table(seqGetData(f, "chromosome"))

# HLA region
seqSetFilterChrom(f, 6, from.bp=29719561, to.bp=32883508)
seqSummary(f, "genotype", check="none")


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
