\name{seqSetFilterChrom}
\alias{seqSetFilterChrom}
\title{Chromosome selection}
\description{
    Select the variants according to the specified chromosome(s).
}
\usage{
seqSetFilterChrom(gdsfile, include=NULL, is.num=NA)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{include}{NULL, or character for specified chromosome(s)}
    \item{is.num}{a logical variable: \code{TRUE}, chromosome code is numeric;
        \code{FALSE}, chromosome is not numeric}
}
\value{
    None.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilter}}
}

\examples{
# the file of GDS
gds.fn <- seqExampleFileName("gds")
# or gds.fn <- "C:/YourFolder/Your_GDS_File.gds"

f <- seqOpen(gds.fn)
seqSummary(f)


seqSetFilterChrom(f, is.num=TRUE)
seqSummary(f, "genotype", check="none")

seqSetFilterChrom(f, is.num=FALSE)
seqSummary(f, "genotype", check="none")

seqSetFilterChrom(f, 1:4)
seqSummary(f, "genotype", check="none")
table(seqGetData(f, "chromosome"))


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
