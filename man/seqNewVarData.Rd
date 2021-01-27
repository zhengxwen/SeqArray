\name{seqNewVarData}
\alias{seqNewVarData}
\title{Variable-length data}
\description{
    Gets a variable-length data object.
}
\usage{
seqNewVarData(len, data)
}
\arguments{
    \item{len}{a non-negative vector for variable lengths}
    \item{data}{a vector of data according to \code{len}}
}
\value{
    Return a \code{SeqVarDataList} object.
}
\details{
    \code{seqGetData()} returns a \code{SeqVarDataList} object for
variable-length data; \code{seqAddValue()} can add a \code{SeqVarDataList}
object to a GDS file.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqGetData}}, \code{\link{seqAddValue}}
}

\examples{
seqNewVarData(c(1,2,1,0,2), c("A", "B", "B", "C", "E", "E"))
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
