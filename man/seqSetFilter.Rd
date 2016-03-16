\name{seqSetFilter-methods}
\alias{seqSetFilter}
\alias{seqSetFilter,SeqVarGDSClass,ANY-method}
\alias{seqSetFilter,SeqVarGDSClass,GRanges-method}
\alias{seqSetFilter,SeqVarGDSClass,GRangesList-method}
\alias{seqSetFilter,SeqVarGDSClass,IRanges-method}
\alias{seqResetFilter}
\title{Set a Filter to Sample or Variant}
\description{
    Sets a filter to sample and/or variant.
}
\usage{
\S4method{seqSetFilter}{SeqVarGDSClass,ANY}(object, variant.sel,
    sample.sel=NULL, variant.id=NULL, sample.id=NULL, samp.sel=NULL,
    action=c("set", "intersect", "push", "push+set", "push+intersect", "pop"),
    verbose=TRUE)
\S4method{seqSetFilter}{SeqVarGDSClass,GRanges}(object, variant.sel,
    rm.txt="chr", verbose=TRUE)
\S4method{seqSetFilter}{SeqVarGDSClass,GRangesList}(object, variant.sel,
    rm.txt="chr", verbose=TRUE)
\S4method{seqSetFilter}{SeqVarGDSClass,IRanges}(object, variant.sel,
    chr, verbose=TRUE)
seqResetFilter(object, sample=TRUE, variant=TRUE, verbose=TRUE)
}
\arguments{
    \item{object}{a \code{\link{SeqVarGDSClass}} object}
    \item{variant.sel}{a logical/raw/index vector indicating the selected
        variants; \code{\link{GRanges}}, a GRanges object for the genomic
        locations; \code{\link{GRangesList}}, a GRangesList object for storing
        a collection of GRanges objects; \code{\link{IRanges}}, a IRanges
        object for storing a collection of range objects}
    \item{sample.sel}{a logical/raw/index vector indicating the selected
        samples}
    \item{variant.id}{ID of selected variants}
    \item{sample.id}{ID of selected samples}
    \item{samp.sel}{it is deprecated since v1.11.12, please use
        \code{sample.sel} instead}
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
    \item{rm.txt}{a character, the characters will be removed from
        \code{seqnames(variant.sel)}}
    \item{chr}{a vector of character for chromsome coding}
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
    \code{\link{seqSetFilterChrom}}, \code{\link{seqGetFilter}},
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
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8)])
set.seed(100)
seqSetFilter(f, variant.id=sample(variant.id, 5))

# get genotypic data
seqGetData(f, "genotype")


## OR
# set sample and variant filters
seqSetFilter(f, sample.sel=c(2,4,6,8))
set.seed(100)
seqSetFilter(f, variant.sel=sample.int(length(variant.id), 5))

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
