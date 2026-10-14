\name{seqSummary}
\alias{seqSummary}
\title{Summarize a SeqArray GDS File}
\description{
    Gets the summary of SeqArray GDS file.
}
\usage{
seqSummary(gdsfile, varname=NULL, check=c("default", "none", "full"),
    verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object, or a file name}
    \item{varname}{if \code{NULL}, check the whole GDS file; or a character
        specifying variable name, and return a description of that variable.
        See \sQuote{Value}}
    \item{check}{should be one of "default", "none", "full";
        \code{check="default"} by default}
    \item{verbose}{if \code{TRUE}, display information}
}
\value{
    If \code{varname=NULL}, a list with components:
    \item{filename}{the file name}
    \item{version}{the version of SeqArray format}
    \item{reference}{genome reference, a character vector (0-length for
        undefined)}
    \item{ploidy}{the number of sets of chromosomes}
    \item{num.sample}{the total number of samples}
    \item{num.variant}{the total number of variants}
    \item{allele}{allele information, see \code{varname="allele"} below}
    \item{annot_qual}{the total number of "annotation/qual" if
        \code{check="none"}, or a summary object including min, max, median,
        mean}
    \item{filter}{filter information, see \code{varname="annotation/filter"}
        below}
    \item{info}{a \code{data.frame} of INFO field: ID, Number, Type,
        Description, Source and Version}
    \item{format}{a \code{data.frame} of FORMAT field: ID, Number, Type and
        Description}
    \item{sample.annot}{a \code{data.frame} of sample annotation with ID,
        Type and Description}

    Otherwise, the return value depends on \code{varname}:
    \describe{
        \item{\code{"genotype"}}{a list with components \code{dim} (an
            integer vector: ploidy, # of samples, # of variants) and
            \code{seldim} (an integer vector: ploidy, # of selected samples,
            # of selected variants); if \code{check="full"}, it also has a
            component \code{digest} for digest validation}
        \item{\code{"chromosome"}}{an integer vector named by chromosome
            codes, giving the number of selected variants on each chromosome
            (numeric codes in ascending order, followed by non-numeric codes,
            e.g., X, Y); if \code{check="full"}, its attribute "digest" is the
            result of digest validation (\code{TRUE}, \code{FALSE}, or
            \code{NA} if no digest). If \code{check="none"}, it returns the
            total number of variants in the file instead, ignoring the variant
            filter}
        \item{\code{"allele"}}{a \code{data.frame} with ID and Description if
            \code{check="none"}; otherwise, a list with components
            \code{value} (that \code{data.frame}) and \code{table} (cross
            tabulation for the number of alleles per site)}
        \item{\code{"$alt"}}{a \code{data.frame} with ID and Description for
            describing the alternative alleles}
        \item{\code{"annotation/filter"}, \code{"$filter"}}{a
            \code{data.frame} with ID and Description if \code{check="none"};
            otherwise, a list with components \code{value} (that
            \code{data.frame}) and \code{table} (cross tabulation for the
            variable 'filter')}
        \item{\code{"annotation/info"}, \code{"$info"}}{a \code{data.frame}
            describing the variables in the folder "annotation/info" with ID,
            Number, Type, Description, Source and Version}
        \item{\code{"annotation/format"}}{a \code{data.frame} describing the
            variables in the folder "annotation/format" with ID, Number, Type
            and Description}
        \item{\code{"$format"}}{a \code{data.frame} describing VCF FORMAT
            header with ID, Number, Type and Description; the first row is
            used for genotypes}
        \item{\code{"sample.annotation"}}{a \code{data.frame} describing
            sample annotation with ID, Type and Description}
        \item{\code{"$reference"}}{the genome reference if it is defined (a
            0-length character vector if undefined)}
        \item{\code{"$contig"}}{the contig information, a \code{data.frame}
            including ID}
        \item{\code{"$digest"}}{a \code{data.frame} with the full names of
            GDS variables, digest codes and validation
            (\code{FALSE}/\code{TRUE})}
    }
}
\details{
    If \code{check="default"}, the function performs regular checking, like
variable dimensions. If \code{check="full"}, it performs more checking, e.g.,
unique sample id, unique variant id, whether genotypic data are in a valid
range or not.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqGetData}}, \code{\link{seqApply}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

seqSummary(gds.fn)

ans <- seqSummary(gds.fn, check="full")
ans

seqSummary(gds.fn, "genotype")
seqSummary(gds.fn, "chromosome")
seqSummary(gds.fn, "allele")
seqSummary(gds.fn, "annotation/filter")
seqSummary(gds.fn, "annotation/info")
seqSummary(gds.fn, "annotation/format")
seqSummary(gds.fn, "sample.annotation")

seqSummary(gds.fn, "$reference")
seqSummary(gds.fn, "$contig")
seqSummary(gds.fn, "$filter")
seqSummary(gds.fn, "$alt")
seqSummary(gds.fn, "$info")
seqSummary(gds.fn, "$format")
seqSummary(gds.fn, "$digest")


# open a GDS file
f <- seqOpen(gds.fn)

# get 'sample.id
samp.id <- seqGetData(f, "sample.id")
# get 'variant.id'
variant.id <- seqGetData(f, "variant.id")

# set sample and variant filters
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8,10)])
set.seed(100)
seqSetFilter(f, variant.id=sample(variant.id, 10))

seqSummary(f, "genotype")
seqSummary(f, "chromosome")

# close a GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
