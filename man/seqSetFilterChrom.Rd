\name{seqSetFilterChrom}
\alias{seqSetFilterChrom}
\title{Chromosome Selection}
\description{
    Selects the variants according to the specified chromosome(s).
}
\usage{
seqSetFilterChrom(gdsfile, include=NULL, is.num=NA, from.bp=NULL, to.bp=NULL,
    verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{include}{NULL, or character for specified chromosome(s)}
    \item{is.num}{a logical variable: \code{TRUE}, chromosome code is numeric;
        \code{FALSE}, chromosome is not numeric; \code{is.num=TRUE} is usually
        used to exclude non-autosomes}
    \item{from.bp}{NULL, no limit; numeric, the lower bound of position}
    \item{to.bp}{NULL, no limit; numeric, the upper bound of position}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    If \code{from.bp} and \code{to.bp} has values, they should be equal-size
as \code{include}. A trio of \code{include}, \code{from.bp} and \code{to.bp}
indicates a region on human genomes. \code{NA} in \code{from.bp} is treated
as 0, and \code{NA} in \code{to.bp} is treated as the maximum of integer
(2^31 - 1).
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

# two regions
seqSetFilterChrom(f, c(1, 6), from.bp=c(1000000, 29719561),
    to.bp=c(90000000, 32883508))
seqSummary(f, "genotype", check="none")
seqGetData(f, "chromosome")

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
