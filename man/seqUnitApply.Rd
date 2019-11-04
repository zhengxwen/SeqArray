\name{seqUnitApply}
\alias{seqUnitApply}
\title{Sliding units of selected variants}
\description{
    Generating units of selected variants via sliding windows.
}
\usage{
seqUnitApply(gdsfile, units, var.name, FUN, as.is=c("none", "list", "unlist"),
    parallel=FALSE, .useraw=FALSE, .progress=FALSE, ...)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{units}{a list of units of selected variants, with class "SeqUnitListClass"}
    \item{var.name}{the variable name(s), see details}
    \item{FUN}{the function to be applied}
    \item{as.is}{returned value: a list, a vector or none;
        \code{as.is} can be a \code{\link{connection}} object, or a GDS node
        \code{\link{gdsn.class}} object; if "unlist" is used, produces a vector
        which contains all the atomic components, via
        \code{unlist(..., recursive=FALSE)}}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details.}
    \item{.useraw}{\code{TRUE}, force to use RAW instead of INTEGER for
        genotypes and dosages; \code{FALSE}, use INTEGER; \code{NA}, use RAW
        instead of INTEGER if possible; for genotypes, 0xFF is missing value
        if RAW is used}
    \item{.progress}{if \code{TRUE}, show progress information}
    \item{...}{optional arguments to \code{FUN}}
}
\value{
    A vector, a list of values or none.
}

\author{Xiuwen Zheng}

\examples{
library(SeqArray)

# open the GDS file
gdsfile <- seqOpen(seqExampleFileName("gds"))

# variant units via sliding windows
units <- seqUnitSlidingWindows(gdsfile)

v1 <- seqUnitApply(gdsfile, units, "genotype", function(x) dim(x)[3L],
    as.is="unlist", .progress=TRUE)
v2 <- seqUnitApply(gdsfile, units, "genotype", function(x) dim(x)[3L],
    as.is="unlist", parallel=2, .progress=TRUE)

all(v1 == lengths(units))
all(v1 == v2)

# close the GDS file
seqClose(gdsfile)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
