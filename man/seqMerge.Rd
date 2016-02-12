\name{seqMerge}
\alias{seqMerge}
\title{Merge Multiple SeqArray GDS Files}
\description{
    Merges multiple SeqArray GDS files.
}
\usage{
seqMerge(gds.fn, out.fn, storage.option="ZIP_RA.default",
    info.var=NULL, fmt.var=NULL, samp.var=NULL, optimize=TRUE, digest=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{gds.fn}{the file names of multiple GDS files}
    \item{out.fn}{the output file name}
    \item{storage.option}{specify the storage and compression options,
        by default \code{\link{seqStorageOption}("ZIP_RA.default")}}
    \item{info.var}{characters, the variable name(s) in the INFO field;
        \code{NULL} for all variables, or \code{character()} excludes all
        INFO variables}
    \item{fmt.var}{characters, the variable name(s) in the FORMAT field;
        \code{NULL} for all variables, or \code{character()} excludes all
        FORMAT variables}
    \item{samp.var}{characters, the variable name(s) in 'sample.annotation';
        or \code{NULL} for all variables}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{digest}{a logical value (TRUE/FALSE) or a character ("md5", "sha1",
        "sha256", "sha384" or "sha512"); add hash codes to the GDS file if TRUE
        or a digest algorithm is specified}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\details{
    The function merges multiple SeqArray GDS files. Users can specify the
compression method and level for the new GDS file. If \code{gds.fn} contains
one file, users can change the storage type to create a new file.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}, \code{\link{seqExport}}
}

\examples{
# the VCF file
vcf.fn <- seqExampleFileName("vcf")

# the number of variants
total.count <- seqVCF_Header(vcf.fn, getnum=TRUE)$num.variant

split.cnt <- 5
start <- integer(split.cnt)
count <- integer(split.cnt)

s <- (total.count+1) / split.cnt
st <- 1L
for (i in 1:split.cnt)
{
    z <- round(s * i)
    start[i] <- st
    count[i] <- z - st
    st <- z
}

fn <- paste0("tmp", 1:split.cnt, ".gds")

# convert to 5 gds files
for (i in 1:split.cnt)
    seqVCF2GDS(vcf.fn, fn[i], start=start[i], count=count[i])

# merge
seqMerge(fn, "tmp.gds")
seqSummary("tmp.gds")


####

vcf.fn <- seqExampleFileName("gds")
file.copy(vcf.fn, "test.gds", overwrite=TRUE)

# modify 'sample.id'
f <- openfn.gds("test.gds", FALSE)
sid <- read.gdsn(index.gdsn(f, "sample.id"))
add.gdsn(f, "sample.id", paste("S", 1:length(sid)), replace=TRUE)
closefn.gds(f)

# merging
seqMerge(c(vcf.fn, "test.gds"), "output.gds")



# delete the temporary files
unlink(c("tmp.gds", "test.gds", "output.gds"), force=TRUE)
unlink(fn, force=TRUE)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
