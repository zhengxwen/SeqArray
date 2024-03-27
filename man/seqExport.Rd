\name{seqExport}
\alias{seqExport}
\title{Export to a GDS File}
\description{
    Exports to a GDS file with selected samples and variants, which are defined
by \code{seqSetFilter()}.
}
\usage{
seqExport(gdsfile, out.fn, info.var=NULL, fmt.var=NULL, samp.var=NULL,
    optimize=TRUE, digest=TRUE, verbose=TRUE, verbose.clean=NA)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{out.fn}{the file name of output GDS file}
    \item{info.var}{characters, the variable name(s) in the INFO field
        for import; or \code{NULL} for all variables}
    \item{fmt.var}{characters, the variable name(s) in the FORMAT field
        for import; or \code{NULL} for all variables}
    \item{samp.var}{characters, the variable name(s) in the folder
        \code{"sample.annotation"}}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{digest}{a logical value (TRUE/FALSE) or a character ("md5", "sha1",
        "sha256", "sha384" or "sha512"); add md5 hash codes to the GDS file
        if TRUE or a digest algorithm is specified}
    \item{verbose}{if \code{TRUE}, show information}
    \item{verbose.clean}{when \code{verbose.clean=NA}, set it to \code{verbose};
        whether display information when calling \code{cleanup.gds} or not;
        only applicable when \code{optimize=TRUE}}
}
\value{
    Return the file name of GDS format with an absolute path.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}, \code{\link{cleanup.gds}}
}

\examples{
# open the GDS file
(gds.fn <- seqExampleFileName("gds"))
(f <- seqOpen(gds.fn))

# get 'sample.id'
head(samp.id <- seqGetData(f, "sample.id"))

# get 'variant.id'
head(variant.id <- seqGetData(f, "variant.id"))

set.seed(100)
# set sample and variant filters
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8,10,12,14,16)])
seqSetFilter(f, variant.id=sample(variant.id, 100))


# export
seqExport(f, "tmp.gds")
seqExport(f, "tmp.gds", info.var=character())
seqExport(f, "tmp.gds", fmt.var=character())
seqExport(f, "tmp.gds", samp.var=character())


# show file
(f1 <- seqOpen("tmp.gds")); seqClose(f1)


# close
seqClose(f)

# delete the temporary file
unlink("tmp.gds")
}

\keyword{gds}
\keyword{VCF}
\keyword{sequencing}
\keyword{genetics}
