\name{seqClose-methods}
\alias{seqClose}
\alias{seqClose,SeqVarGDSClass-method}
\alias{seqClose,gds.class-method}
\title{Close the SeqArray GDS File}
\description{
    Closes a sequence GDS file which is open.
}

\usage{
\S4method{seqClose}{gds.class}(object)
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
        \item{\code{\link{gds.class}}, close a general GDS file}
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
