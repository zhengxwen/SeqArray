\name{seqSystem}
\alias{seqSystem}
\title{Get the parameters in the GDS system}
\description{
    Get a list of parameters in the GDS system
}

\usage{
seqSystem()
}
\value{
    A list including
    \item{num.logical.core}{the number of logical cores}
    \item{compiler.flag}{SIMD instructions supported by the compiler}
    \item{options}{list all options associated with SeqArray GDS format or
        packages}
}

\author{Xiuwen Zheng}

\examples{
seqSystem()
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
