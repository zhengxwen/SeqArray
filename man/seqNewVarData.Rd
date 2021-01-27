\name{seqNewVarData}
\alias{seqNewVarData}
\alias{seqListVarData}
\title{Variable-length data}
\description{
    Gets a variable-length data object.
}
\usage{
seqNewVarData(len, data)
seqListVarData(obj)
}
\arguments{
    \item{len}{a non-negative vector for variable lengths}
    \item{data}{a vector of data according to \code{len}}
    \item{obj}{a \code{SeqVarDataList} object}
}
\value{
    Return a \code{SeqVarDataList} object.
}
\details{
    \code{seqNewVarData()} creates a \code{SeqVarDataList} object for variable-
length data, and \code{seqListVarData()} converts the \code{SeqVarDataList}
object to a list.
    \code{seqGetData()} returns a \code{SeqVarDataList} object for
variable-length data; \code{seqAddValue()} can add a \code{SeqVarDataList}
object to a GDS file.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqGetData}}, \code{\link{seqAddValue}}
}

\examples{
obj <- seqNewVarData(c(1,2,1,0,2), c("A", "B", "B", "C", "E", "E"))
obj

seqListVarData(obj)
}

\keyword{gds}
\keyword{sequencing}
