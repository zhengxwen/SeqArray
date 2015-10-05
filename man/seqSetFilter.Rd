\name{seqSetFilter-methods}
\alias{seqSetFilter}
\alias{seqSetFilter,SeqVarGDSClass-method}
\alias{seqResetFilter}
\title{Set a Filter to Sample or Variant}
\description{
    Sets a filter to sample and/or variant.
}
\usage{
\S4method{seqSetFilter}{SeqVarGDSClass}(object,
    sample.id=NULL, variant.id=NULL, samp.sel=NULL, variant.sel=NULL,
    action=c("set", "intersect", "push", "push+set", "push+intersect", "pop"),
    verbose=TRUE)

seqResetFilter(object, sample=TRUE, variant=TRUE, verbose=TRUE)
}
\arguments{
    \item{object}{a \code{\link{SeqVarGDSClass}} object}
    \item{sample.id}{IDs of selected samples}
    \item{variant.id}{IDs of selected variants}
    \item{samp.sel}{a logical/raw vector indicating the selected samples}
    \item{variant.sel}{a logical/raw vector indicating the selected variants}
    \item{action}{\code{"set"} -- set the current filter via \code{sample.id},
            \code{variant.id}, \code{samp.sel} or \code{variant.sel};
        \code{"intersect"} -- set the current filter to the intersection of
            selected samples and/or variants;
        \code{"push"} -- push the current filter to the stack, and it could
            be recovered by \code{"pop"} later, no change on the current
            filter;
        \code{"push+set"} -- push the current filter to the stack, and changes
            the current filter via \code{sample.id}, \code{variant.id},
            \code{samp.sel} or \code{variant.sel};
        \code{"push+intersect"} -- push the current filter to the stack, and
            set the current filter to the intersection of selected samples
            and/or variants;
        \code{"pop"} -- pop up the last filter}
    \item{sample}{logical, if \code{TRUE}, include all samples}
    \item{variant}{logical, if \code{TRUE}, include all variants}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    \code{seqResetFilter(file)} is equivalent to \code{seqSetFilter(file)},
where the selection arguments in \code{seqSetFilter} are \code{NULL}.
}
\value{
    None.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqGetFilter}}, \code{\link{seqSetFilterChrom}},
    \code{\link{seqGetData}}, \code{\link{seqApply}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
(f <- seqOpen(gds.fn))

# get 'sample.id
(samp.id <- seqGetData(f, "sample.id"))
# "NA06984" "NA06985" "NA06986" ...

# get 'variant.id'
head(variant.id <- seqGetData(f, "variant.id"))

# get 'chromosome'
table(seqGetData(f, "chromosome"))

# get 'allele'
head(seqGetData(f, "allele"))
# "T,C" "G,A" "G,A" ...


# set sample and variant filters
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8,10)])
set.seed(100)
seqSetFilter(f, variant.id=sample(variant.id, 10))

# get genotypic data
seqGetData(f, "genotype")



## set the intersection

seqResetFilter(f)
seqSetFilterChrom(f, 10L)
seqSummary(f, "genotype", check="none")

AF <- seqAlleleFreq(f)
table(AF <= 0.9)

seqSetFilter(f, variant.sel=(AF<=0.9), action="intersect")
seqSummary(f, "genotype", check="none")


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
