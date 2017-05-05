\name{seqAsVCF}
\alias{seqAsVCF}
\alias{seqAsVCFInit}
\title{VariantAnnotation objects}
\description{
    Imports the VariantAnnotation package and export the generic methods
associated with SeqVarGDSClass class.
}
\usage{
seqAsVCFInit()
seqAsVCF(gdsfile, ...)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{...}{optional arguments to \code{asVCF}}
}
\details{
    The generic methods \code{ref}, \code{alt}, \code{qual}, \code{filt},
\code{fixed}, \code{header}, \code{info}, \code{geno}, \code{asVCF},
\code{rowRanges} and \code{colData} are exported and associated with a
\code{\link{SeqVarGDSClass}} object after calling \code{seqAsVCFInit()}
or \code{seqAsVCF()}.
}
\value{
    None or a \code{\link[VariantAnnotation]{CollapsedVCF}} object.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{asVCF,SeqVarGDSClass-method}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))

\dontrun{
seqAsVCFInit()  ## enable exporting the generic methods
seqAsVCF(gds)
}

seqClose(gds)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
