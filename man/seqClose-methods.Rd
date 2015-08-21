\name{seqClose-methods}
\alias{seqClose}
\alias{seqClose,SeqVarGDSClass-method}

\title{Close the SeqArray GDS File}
\description{
    Closes a sequence GDS file.
}

\usage{
\S4method{seqClose}{SeqVarGDSClass}(object)
}
\arguments{
    \item{object}{a SeqArray object}
}
\value{
    None.
}
\details{
    If \code{object} is
    \itemize{
        \item{\code{\link{SeqVarGDSClass}}, close the sequence GDS file.}
    }
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqOpen}}
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
