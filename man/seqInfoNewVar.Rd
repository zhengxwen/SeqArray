\name{seqInfoNewVar}
\alias{seqInfoNewVar}
\title{Add a variable to the INFO field}
\description{
	Add a new variable to the INFO field in the specified GDS file.
}
\usage{
seqInfoNewVar(gdsfile, var.name, variant.id, val,
    description="", compress=c("ZIP.MAX", ""), no.data.index=TRUE)
}
\arguments{
	\item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
	\item{var.name}{the variable name(s), see details}
	\item{variant.id}{the variant ID(s), should have the same order as
		IDs in the GDS file}
	\item{val}{a vector or a matrix with the same order as \code{variant.id};
		if it is a matrix, the number of columns should be equal to the
		length of \code{variant.id}, see details}
	\item{description}{the variable description}
	\item{compress}{to specify the compression algorithm: \code{""}
		for no compression; see \code{compress} in the function
		\code{\link[gdsfmt]{add.gdsn}}}
	\item{no.data.index}{applicable only if \code{val} is a numeric vector
		or a factor variable; if \code{no.data.index=TRUE}, non-exist values
		in \code{val} will be replaced by \code{NA} or \code{NaN} with respect
		to the variants not specified in \code{variant.id}, and no index
		data associated with this variable are created}
}
\details{
	The variable name should be \code{"sample.id"}, \code{"variant.id"},
\code{"position"}, \code{"chromosome"}, \code{"allele"}, \code{"annotation/id"},
\code{"annotation/qual"}, \code{"annotation/filter"},
\code{"annotation/info/VARIABLE_NAME"}, or \code{"annotation/format/VARIABLE_NAME"}.

	The argument \code{val} should be integers, numeric values, a logical
variable, characters or factors. If \code{val} is a logical variable, one-bit
storage mode will be used to store this variable, which corresponds to the
variable defined with 'Type=Flag' in the VCF format.
}
\value{
	None.
}

\author{Xiuwen Zheng}
\seealso{
	\code{\link{seqGetData}}
}

\examples{
# the file of GDS
gds.fn <- seqExampleFileName("gds")

file.copy(gds.fn, "test.gds", overwrite=TRUE)

# display
(f <- seqOpen("test.gds", readonly=FALSE))

# get variant IDs
variant.id <- seqGetData(f, "variant.id")


####    add variables to the INFO field    ####

set.seed(100)

seqInfoNewVar(f, "int", variant.id[1:15], sample.int(256, 15, TRUE),
	"integer variable")
seqGetData(f, "annotation/info/int")

seqInfoNewVar(f, "int.2", variant.id[1:15], sample.int(256, 15, TRUE),
	"integer variable", no.data.index=FALSE)
seqGetData(f, "annotation/info/int.2")


seqInfoNewVar(f, "numeric", variant.id[15:30], rnorm(16),
	"numeric variable")
seqGetData(f, "annotation/info/numeric")

seqInfoNewVar(f, "numeric.2", variant.id[15:30], rnorm(16),
	"numeric variable", no.data.index=FALSE)
seqGetData(f, "annotation/info/numeric.2")


seqInfoNewVar(f, "flag", variant.id[4:9], rep(c(FALSE, TRUE), 3),
	"flag variable")
# stored in `bit1'
seqGetData(f, "annotation/info/flag")


seqInfoNewVar(f, "factor", variant.id,
	factor(c("ABC", "DDD", "CVX")[sample(1:3, length(variant.id), TRUE)]),
	"string/factor variable")
# stored in `int32' with attributes
seqGetData(f, "annotation/info/factor")


seqInfoNewVar(f, "string", variant.id,
	c("ABC", "DDD", "CVX")[sample(1:3, length(variant.id), TRUE)],
	"string variable")
seqGetData(f, "annotation/info/string")


# show the file
f


# the corresponding VCF file
seqGDS2VCF(f, "test.vcf.gz")
txt <- strsplit(readLines("test.vcf.gz", n=40), "\t")[-c(1:21)]

# the INFO field:
sapply(txt, function(x) x[8])


# close the GDS file
seqClose(f)


# delete the temporary files
unlink("test.gds", force=TRUE)
unlink("test.vcf.gz", force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
