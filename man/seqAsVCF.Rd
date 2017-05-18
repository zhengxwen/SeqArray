\name{seqAsVCF}
\alias{seqAsVCF}
\title{VariantAnnotation objects}
\description{
    Create a \link[VariantAnnotation]{VCF-class} object
}
\usage{
seqAsVCF(x, chr.prefix="", info=NULL, geno=NULL)
}
\arguments{
    \item{x}{a \code{\link{SeqVarGDSClass}} object}
    \item{chr.prefix}{prefix to add to seqlevels}
    \item{info}{which INFO fields to return}
    \item{geno}{which GENO fields to return}
}
\details{
Coerces a SeqVarGDSClass object to a
\link[VariantAnnotation]{VCF-class} object.
Row names correspond to the variant.id.  \code{info} and
\code{geno} specify the 'INFO' and 'GENO' (FORMAT) fields to
return, respectively.  If not specified, all fields are returned;
if 'NA' no fields are returned.  Use \code{\link{seqSetFilter}}
prior to calling \code{seqAsVCF} to specify samples and variants to
return.

The \pkg{VariantAnnotation} package should be loaded to explore this object.
}
\value{
    A \code{\link[VariantAnnotation]{CollapsedVCF}} object.
}

\author{Stephanie Gogarten, Xiuwen Zheng}
\seealso{
    \code{\link[VariantAnnotation]{VCF}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))

\dontrun{
library(VariantAnnotation)
seqAsVCF(gds)
}

seqClose(gds)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
