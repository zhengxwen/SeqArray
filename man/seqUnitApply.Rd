\name{seqUnitApply}
\alias{seqUnitApply}
\title{Sliding units of selected variants}
\description{
    Generating units of selected variants via sliding windows.
}
\usage{
seqUnitApply(gdsfile, units, var.name, FUN, as.is=c("none", "list", "unlist"),
    parallel=FALSE, ..., .bl_size=256L, .progress=FALSE, .useraw=FALSE, .padNA=TRUE,
    .envir=NULL)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{units}{a list of units of selected variants, with S3 class "SeqUnitListClass"}
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
    \item{.bl_size}{chuck size, the increment for load balancing, 256 for units}
    \item{.progress}{if \code{TRUE}, show progress information}
    \item{.useraw}{\code{TRUE}, force to use RAW instead of INTEGER for
        genotypes and dosages; \code{FALSE}, use INTEGER; \code{NA}, use RAW
        instead of INTEGER if possible; for genotypes, 0xFF is missing value
        if RAW is used}
    \item{.padNA}{\code{TRUE}, pad a variable-length vector with NA if the number of data
        points for each variant is not greater than 1}
    \item{.envir}{NULL, an environment object, or a list/data.frame}
    \item{...}{optional arguments to \code{FUN}}
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

\code{"$num_allele"} returns an integer vector with the numbers of distinct
alleles.

\code{"$ref"} returns a character vector of reference alleles

\code{"$alt"} returns a character vector of alternative alleles (delimited by
comma)

\code{"$chrom_pos"} returns characters with the combination of chromosome and
position, e.g., "1:1272721". \code{"$chrom_pos_allele"} returns characters with
the combination of chromosome, position and alleles, e.g., "1:1272721_A_G"
(i.e., chr:position_REF_ALT).

\code{"$variant_index"} returns the indices of selected variants starting from 1.
\code{"$sample_index"} returns the indices of selected samples starting from 1.
}
\value{
    A vector, a list of values or none.
}

\author{Xiuwen Zheng}

\examples{
# open the GDS file
gdsfile <- seqOpen(seqExampleFileName("gds"))

# variant units via sliding windows
units <- seqUnitSlidingWindows(gdsfile)

v1 <- seqUnitApply(gdsfile, units, "genotype", function(x) dim(x)[3L],
    as.is="unlist", .progress=TRUE)
v2 <- seqUnitApply(gdsfile, units, "genotype", function(x) dim(x)[3L],
    as.is="unlist", parallel=2, .progress=TRUE)

all(v1 == lengths(units$index))
all(v1 == v2)


# call with an external R variable
ext <- list(x=1:1348/10)
v3 <- seqUnitApply(gdsfile, units, "$:x", function(x) x,
    as.is="list", .progress=TRUE, .envir=ext)
head(units$index)
head(v3)

table(sapply(seq_along(units$index), function(i) all(units$index[[i]] == v3[[i]]*10)))
# all TRUE


# close the GDS file
seqClose(gdsfile)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
