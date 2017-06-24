\name{seqGetData}
\alias{seqGetData}
\title{Get Data}
\description{
    Gets data from a SeqArray GDS file.
}
\usage{
seqGetData(gdsfile, var.name, .useraw=FALSE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{var.name}{the variable name, see details}
    \item{.useraw}{\code{TRUE}, force to use RAW instead of INTEGER for
        genotypes and dosages; \code{FALSE}, use INTEGER; \code{NA}, use RAW
        for small numbers instead of INTEGER if possible; 0xFF is missing value
        if RAW is used}
}
\value{
    Return vectors or lists.
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

\code{"$chrom_pos"} returns characters with the combination of chromosome and
position, e.g., "1_1272721".

\code{"$dosage"} returns a RAW/INTEGER matrix for the dosages of reference
allele.

\code{"$num_allele"} returns an integer vector with the numbers of distinct
alleles.

\code{"$ref"} returns a character vector of reference alleles

\code{"$alt"} returns a character vector of alternative alleles (delimited by
comma)
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqSetFilter}}, \code{\link{seqApply}}
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

# get '$num_allele'
head(seqGetData(f, "$num_allele"))


# set sample and variant filters
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8,10)])
set.seed(100)
seqSetFilter(f, variant.id=sample(variant.id, 10))

# get genotypic data
seqGetData(f, "genotype")

# get annotation/info/DP
seqGetData(f, "annotation/info/DP")

# get annotation/info/AA, a variable-length dataset
seqGetData(f, "annotation/info/AA")
# $length              <- indicating the length of each variable-length data
# [1] 1 1 1 1 1 1 ...
# $data                <- the data according to $length
# [1] "T" "C" "T" "C" "G" "C" ...

# get annotation/format/DP, a variable-length dataset
seqGetData(f, "annotation/format/DP")
# $length              <- indicating the length of each variable-length data
# [1] 1 1 1 1 1 1 ...
# $data                <- the data according to $length
#      variant
# sample [,1] [,2] [,3] [,4] [,5] [,6] ...
#  [1,]   25   25   22    3    4   17  ...

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
