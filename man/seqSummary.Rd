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
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{varname}{if \code{NULL}, check the whole GDS file; or a character
        specifying variable name, and return a description of that variable.
        See details}
    \item{check}{should be one of "default", "none", "full"}
    \item{verbose}{if \code{TRUE}, display information}
}
\value{
    If \code{varname=NULL}, the function returns a list:
    \item{filename}{the file name}
    \item{version}{the version of SeqArray format}
    \item{reference}{genome reference, a character vector (0-length for
        undefined)}
    \item{ploidy}{the number of sets of chromosomes}
    \item{num.sample}{the total number of samples}
    \item{num.variant}{the total number of variants}
    \item{allele}{allele information, see \code{seqSummary(gdsfile, "allele")}}
    \item{annot_qual}{the total number of "annotation/qual" if
        \code{check="none"}, or a summary object including min, max, median,
        mean}
    \item{filter}{filter information, see
        \code{seqSummary(gdsfile, "annotation/filter")}}
    \item{info}{a \code{data.frame} of INFO field: ID, Number, Type,
        Description, Source and Version}
    \item{format}{a \code{data.frame} of FORMAT field: ID, Number, Type and
        Description}
    \item{sample.annot}{a \code{data.frame} of sample annotation with ID,
        Type and Description}

    ---
    \code{seqSummary(gdsfile, "genotype", check="none", verbose=FALSE)} returns
a list with components:
    \item{dim}{ploidy, # of samples, # of variants}
    \item{seldim}{# of selected samples, # of selected variants}

    ---
    \code{seqSummary(gdsfile, "allele")} returns a data.frame with ID and
descriptions (\code{check="none"}), or a list with components:
    \item{value}{a data.frame with ID and Description}
    \item{table}{cross tabulation for the number of alleles per site}

    ---
    \code{seqSummary(gdsfile, "$alt")} returns a data.frame with ID and
Description for describing the alternative alleles.

    ---
    \code{seqSummary(gdsfile, "annotation/filter")} or
\code{seqSummary(gdsfile, "$filter")} returns a data.frame with ID and
description (\code{check="none"}), or a list with components:
    \item{value}{a data.frame with ID and Description}
    \item{table}{cross tabulation for the variable 'filter'}

    ---
    \code{seqSummary(gdsfile, "annotation/info")} or
\code{seqSummary(gdsfile, "$info")} returns a \code{data.frame} describing
the variables in the folder "annotation/info" with ID, Number, Type,
Description, Source and Version.

    ---
    \code{seqSummary(gdsfile, "annotation/format")} returns a \code{data.frame}
describing the variables in the folder "annotation/format" with ID, Number,
Type and Description.

    ---
    \code{seqSummary(gdsfile, "sample.annotation")} returns a \code{data.frame}
describing sample annotation with ID, Type and Description.

    ---
    \code{seqSummary(gdsfile, "$reference")} returns the genome reference
if it is defined (a 0-length character vector if undefined).

    ---
    \code{seqSummary(gdsfile, "$contig")} returns the contig information,
a \code{data.frame} including ID.

    ---
    \code{seqSummary(gdsfile, "$format")} returns a \code{data.frame}
describing VCF FORMAT header with ID, Number, Type and Description. The first
row is used for genotypes.

    ---
    \code{seqSummary(gdsfile, "$digest")} returns a \code{data.frame} with
the full names of GDS variables, digest codes and validation (FALSE/TRUE).
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

# close a GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
