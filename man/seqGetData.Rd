\name{seqGetData}
\alias{seqGetData}
\title{Get Data}
\description{
    Gets data from a SeqArray GDS file.
}
\usage{
seqGetData(gdsfile, var.name, .useraw=FALSE, .padNA=TRUE, .tolist=FALSE,
    .envir=NULL)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{var.name}{a variable name or a character vector, see details;
        if \code{character()}, return \code{NULL}}
    \item{.useraw}{\code{TRUE}, force to use RAW instead of INTEGER for
        genotypes and dosages; \code{FALSE}, use INTEGER; \code{NA}, use RAW
        for small numbers instead of INTEGER if possible; 0xFF is missing value
        if RAW is used}
    \item{.padNA}{\code{TRUE}, pad a variable-length vector with NA if the
        number of data points for each variant is not greater than 1}
    \item{.tolist}{if \code{TRUE}, return a list of vectors instead of the
        structure \code{list(length, data)} for variable-length data}
    \item{.envir}{\code{NULL}, an environment object, a list or a
        \code{data.frame}}
}
\value{
    Return vectors, matrices or lists (with \code{length} and \code{data}
components) with a class name \code{SeqVarDataList}.
}
\details{
    The variable name should be \code{"sample.id"}, \code{"variant.id"},
\code{"position"}, \code{"chromosome"}, \code{"allele"}, \code{"genotype"},
\code{"annotation/id"}, \code{"annotation/qual"}, \code{"annotation/filter"},
\code{"annotation/info/VARIABLE_NAME"}, or
\code{"annotation/format/VARIABLE_NAME"}.

\code{"@genotype"}, \code{"annotation/info/@VARIABLE_NAME"} or
\code{"annotation/format/@VARIABLE_NAME"} are used to obtain the index
associated with these variables.

\code{"$dosage"} is also allowed for the dosages of reference allele (integer:
0, 1, 2 and NA for diploid genotypes).

\code{"$dosage_alt"} returns a RAW/INTEGER matrix for the dosages of alternative
allele without distinguishing different alternative alleles.
\code{"$dosage_alt2"} allow the alleles are partially missing (e.g., genotypes
on chromosome X for males)

\code{"$dosage_sp"} returns a sparse matrix (dgCMatrix) for the dosages of
alternative allele without distinguishing different alternative alleles.
\code{"$dosage_sp2"} allow the alleles are partially missing (e.g., genotypes
on chromosome X for males)

\code{"$num_allele"} returns an integer vector with the numbers of distinct
alleles.

\code{"$ref"} returns a character vector of reference alleles.
\code{"$alt"} returns a character vector of alternative alleles (delimited by
comma).

\code{"$chrom_pos"} returns characters with the combination of chromosome and
position, e.g., "1:1272721". \code{"$chrom_pos2"} is similar to
\code{"$chrom_pos"}, except the suffix "_1" is added to the first duplicate
following the variant, "_2" is added to the second duplicate, and so on.
\code{"$chrom_pos_allele"} returns characters with the combination of
chromosome, position and alleles, e.g., "1:1272721_A_G"
(i.e., chr:position_REF_ALT).

\code{"$variant_index"} returns the indices of selected variants starting
from 1, and \code{"$sample_index"} returns the indices of selected samples
starting from 1.

\code{"$:VAR"} return the variable "VAR" from \code{.envir} according to the
selected variants.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilter}}, \code{\link{seqApply}},
    \code{\link{seqNewVarData}}, \code{\link{seqListVarData}}
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

# get '$chrom_pos'
head(seqGetData(f, "$chrom_pos"))

# get '$dosage'
seqGetData(f, "$dosage")[1:6, 1:10]

# get a sparse matrix of dosages
seqGetData(f, "$dosage_sp")[1:6, 1:10]

# get '$num_allele'
head(seqGetData(f, "$num_allele"))


# set sample and variant filters
set.seed(100)
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8,10)])
seqSetFilter(f, variant.id=sample(variant.id, 10))

# get a list
seqGetData(f, c(chr="chromosome", pos="position", allele="allele"))

# get the indices of selected variants/samples
seqGetData(f, "$variant_index")
seqGetData(f, "$sample_index")

# get genotypic data
seqGetData(f, "genotype")

# get annotation/info/DP
seqGetData(f, "annotation/info/DP")

# get annotation/info/AA, a variable-length dataset
seqGetData(f, "annotation/info/AA", .padNA=FALSE)
# $length              <- indicating the length of each variable-length data
# [1] 1 1 1 1 1 1 ...
# $data                <- the data according to $length
# [1] "T" "C" "T" "C" "G" "C" ...

# or return a simplified vector
seqGetData(f, "annotation/info/AA", .padNA=TRUE)

# get annotation/format/DP, a variable-length dataset
seqGetData(f, "annotation/format/DP")
# $length              <- indicating the length of each variable-length data
# [1] 1 1 1 1 1 1 ...
# $data                <- the data according to $length
#      variant
# sample [,1] [,2] [,3] [,4] [,5] [,6] ...
#  [1,]   25   25   22    3    4   17  ...


# get values from R environment
env <- new.env()
env$x <- 1:1348 / 10
env$x[seqGetData(f, "$variant_index")]
seqGetData(f, "$:x", .envir=env)


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
