\name{seqSetFilter-methods}
\alias{seqSetFilter}
\alias{seqSetFilter,SeqVarGDSClass,ANY-method}
\alias{seqSetFilter,SeqVarGDSClass,GRanges-method}
\alias{seqSetFilter,SeqVarGDSClass,GRangesList-method}
\alias{seqSetFilter,SeqVarGDSClass,IRanges-method}
\alias{seqResetFilter}
\alias{seqSetFilterChrom}
\alias{seqSetFilterPos}
\title{Set a Filter to Sample or Variant}
\description{
    Sets a filter to sample and/or variant.
}
\usage{
\S4method{seqSetFilter}{SeqVarGDSClass,ANY}(object, variant.sel,
    sample.sel=NULL, variant.id=NULL, sample.id=NULL,
    action=c("set", "intersect", "push", "push+set", "push+intersect", "pop"),
    verbose=TRUE)
\S4method{seqSetFilter}{SeqVarGDSClass,GRanges}(object, variant.sel,
    rm.txt="chr", intersect=FALSE, verbose=TRUE)
\S4method{seqSetFilter}{SeqVarGDSClass,GRangesList}(object, variant.sel,
    rm.txt="chr", intersect=FALSE, verbose=TRUE)
\S4method{seqSetFilter}{SeqVarGDSClass,IRanges}(object, variant.sel,
    chr, intersect=FALSE, verbose=TRUE)
seqResetFilter(object, sample=TRUE, variant=TRUE, verbose=TRUE)
seqSetFilterChrom(object, include=NULL, is.num=NA, from.bp=NULL, to.bp=NULL,
    intersect=FALSE, verbose=TRUE)
seqSetFilterPos(object, chr, pos, intersect=FALSE, multi.pos=FALSE,
    verbose=TRUE)
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
    \item{pos}{a vector of numeric values for genome coordinate}
    \item{sample}{logical, if \code{TRUE}, include all samples}
    \item{variant}{logical, if \code{TRUE}, include all variants}
    \item{include}{NULL, or a vector of characters for specified chromosome(s)}
    \item{is.num}{a logical variable: \code{TRUE}, chromosome code is numeric;
        \code{FALSE}, chromosome is not numeric; \code{is.num=TRUE} is usually
        used to exclude non-autosomes}
    \item{from.bp}{NULL, no limit; a numeric vector, the lower bound of
    	position}
    \item{to.bp}{NULL, no limit; a numeric vector, the upper bound of
    	position}
    \item{intersect}{if \code{FALSE}, the candidate samples/variants for
        selection are all samples/variants (by default); if \code{TRUE}, the
        candidate samples/variants are from the selected samples/variants
        defined via the previous call}
    \item{multi.pos}{\code{FALSE}, use the first matched position;
        \code{TRUE}, allow multiple variants at the same position}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    \code{seqResetFilter(file)} is equivalent to \code{seqSetFilter(file)},
where the selection arguments in \code{seqSetFilter} are \code{NULL}.

    If \code{from.bp} and \code{to.bp} has values, they should be equal-size
as \code{include}. A trio of \code{include}, \code{from.bp} and \code{to.bp}
indicates a region on human genomes. \code{NA} in \code{from.bp} is treated
as 0, and \code{NA} in \code{to.bp} is treated as the maximum of integer
(2^31 - 1).
}
\value{
    None.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilterChrom}}, \code{\link{seqSetFilterCond}},
    \code{\link{seqGetFilter}}, \code{\link{seqGetData}},
    \code{\link{seqApply}}
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



## chromosome

seqResetFilter(f)

seqSetFilterChrom(f, is.num=TRUE)
seqSummary(f, "genotype", check="none")

seqSetFilterChrom(f, is.num=FALSE)
seqSummary(f, "genotype", check="none")

seqSetFilterChrom(f, 1:4)
seqSummary(f, "genotype", check="none")
table(seqGetData(f, "chromosome"))

# HLA region
seqSetFilterChrom(f, 6, from.bp=29719561, to.bp=32883508)
seqSummary(f, "genotype", check="none")

# two regions
seqSetFilterChrom(f, c(1, 6), from.bp=c(1000000, 29719561),
    to.bp=c(90000000, 32883508))
seqSummary(f, "genotype", check="none")
seqGetData(f, "chromosome")


## intersection option

seqResetFilter(f)
seqSetFilterChrom(f, 6, from.bp=29719561, to.bp=32883508)  # MHC
seqSetFilterChrom(f, include=6)  # chromosome 6

seqResetFilter(f)
seqSetFilterChrom(f, 6, from.bp=29719561, to.bp=32883508)  # MHC
seqSetFilterChrom(f, include=6, intersect=TRUE)  # MHC region only



# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
