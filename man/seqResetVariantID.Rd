\name{seqResetVariantID}
\alias{seqResetVariantID}
\title{Reset Variant ID in SeqArray GDS Files}
\description{
    Resets the variant IDs in multiple SeqArray GDS files.
}
\usage{
seqResetVariantID(gds.fn, digest=TRUE, optimize=TRUE, verbose=TRUE)
}
\arguments{
    \item{gds.fn}{a list of file names of multiple GDS files}
    \item{digest}{a logical value, if TRUE, add a md5 hash code}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    None.
}
\details{
    The variant IDs will be replaced by the numbers in sequential order and adjacent to
each file. The variant ID starts from 1 in the first GDS file.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
fn <- seqExampleFileName("gds")

file.copy(fn, "tmp1.gds", overwrite=TRUE)
file.copy(fn, "tmp2.gds", overwrite=TRUE)

gds.fn <- c("tmp1.gds", "tmp2.gds")
seqResetVariantID(gds.fn)

f <- seqOpen("tmp1.gds")
head(seqGetData(f, "variant.id"))
seqClose(f)

f <- seqOpen("tmp2.gds")
head(seqGetData(f, "variant.id"))
seqClose(f)


# delete the temporary files
unlink(gds.fn, force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
