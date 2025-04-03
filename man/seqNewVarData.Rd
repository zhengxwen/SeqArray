\name{seqNewVarData}
\alias{seqNewVarData}
\alias{seqListVarData}
\title{Variable-length data}
\description{
    Gets a variable-length data object.
}
\usage{
seqNewVarData(len, data)
seqListVarData(obj, useList=FALSE)
}
\arguments{
    \item{len}{a non-negative vector for variable lengths}
    \item{data}{a vector of data according to \code{len}}
    \item{obj}{a \code{SeqVarDataList} object or a compressed list (defined
        in \pkg{IRanges})}
    \item{useList}{if \code{TRUE}, return a compressed List defined in
        \pkg{IRanges} (e.g., \code{\link{seqGetData}})}
}
\value{
    Return a \code{SeqVarDataList} object or a \code{CompressedAtomicList}
object (e.g., \code{IntegerList()}).
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
    \code{\link{seqGetData}}, \code{\link{seqAddValue}},
    \code{\link{IntegerList}}
}

\examples{
obj <- seqNewVarData(c(1,2,1,0,2), c("A", "B", "B", "C", "E", "E"))
obj

seqListVarData(obj)  # a list
(a <- seqListVarData(obj, useList=TRUE))  # CharacterList
seqListVarData(a)  # a list
}

\keyword{gds}
\keyword{sequencing}
