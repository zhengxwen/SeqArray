\name{seqVCF2GDS}
\alias{seqVCF2GDS}
\title{Reformat VCF Files}
\description{
    Reformats Variant Call Format (VCF) files.
}
\usage{
seqVCF2GDS(vcf.fn, out.fn, header=NULL, genotype.var.name="GT",
    genotype.storage=c("bit2", "bit4", "bit8"),
    storage.option=seqStorage.Option(),
    info.import=NULL, fmt.import=NULL, ignore.chr.prefix="chr",
    reference=NULL, start=1L, count=-1L, optimize=TRUE, raise.error=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{vcf.fn}{the file name(s) of VCF format}
    \item{out.fn}{the file name of output GDS file}
    \item{header}{if NULL, \code{header} is set to be
        \code{\link{seqVCF.Header}(vcf.fn)}}
    \item{genotype.var.name}{the ID for genotypic data in the FORMAT column;
        "GT" by default, VCFv4.0}
    \item{genotype.storage}{"bit2" by default; with respect to the compression
        size and access speed, "bit2" is the most efficient when most of
        variants are biallelic.}
    \item{storage.option}{specify the storage and compression options,
        by default \code{\link{seqStorage.Option}}, see details}
    \item{info.import}{characters, the variable name(s) in the INFO field
        for import; or \code{NULL} for all variables}
    \item{fmt.import}{characters, the variable name(s) in the FORMAT field
        for import; or \code{NULL} for all variables}
    \item{ignore.chr.prefix}{a vector of character, indicating the prefix of
        chromosome which should be ignored, like "chr"; it is not
        case-sensitive}
    \item{reference}{genome reference, like "hg19", "GRCh37"; if the genome
        reference is not available in VCF files, users could specify the
        reference here}
    \item{start}{the starting variant if importing part of VCF files}
    \item{count}{the maximum count of variant if importing part of VCF files,
        -1 indicates importing to the end}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{raise.error}{\code{TRUE}: throw an error if numeric conversion fails;
        \code{FALSE}: get missing value if numeric conversion fails}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\details{
    GDS -- Genomic Data Structures used for storing genetic array-oriented
        data, and the file format defined in the \link{gdsfmt} package.

    VCF -- The Variant Call Format (VCF), which is a generic format for
        storing DNA polymorphism data such as SNPs, insertions, deletions and
        structural variants, together with rich annotations.

    If there are more than one files in \code{vcf.fn}, \code{seqVCF2GDS} will
merge all VCF files together if they contain the same samples. It is useful
to merge genomic variants if VCF data are divided by chromosomes.

    The real numbers in the VCF file(s) are stored in 32-bit floating-point
format by default. Users can set \code{seqStorage.Option(float.mode="float64")}
to switch to 64-bit floating point format. Or packed real numbers can be
adopted by setting
\code{seqStorage.Option(float.mode="packedreal16:scale=0.0001")}.

    By default, the compression method is "ZIP_RA" (zlib algorithm with default
compression level + multiple independent data blocks). Users can maximize the
compression ratio by \code{seqStorage.Option("ZIP_RA.max")}.
    LZ4 (an extremely fast compression algorithm,
\url{http://cyan4973.github.io/lz4/}) is an option via
\code{storage.option=seqStorage.Option("LZ4_RA")}.
}
\references{
    The variant call format and VCFtools.
    Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA,
    Handsaker RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R;
    1000 Genomes Project Analysis Group.
    Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF.Header}}, \code{\link{seqStorage.Option}},
    \code{\link{seqMerge}}, \code{\link{seqGDS2VCF}}
}

\examples{
# the VCF file
vcf.fn <- seqExampleFileName("vcf")

# convert
seqVCF2GDS(vcf.fn, "tmp.gds")

# display
(f <- seqOpen("tmp.gds"))
seqClose(f)



# convert without the INFO fields
seqVCF2GDS(vcf.fn, "tmp.gds", info.import=character(0))

# display
(f <- seqOpen("tmp.gds"))
seqClose(f)



# convert without the INFO and FORMAT fields
seqVCF2GDS(vcf.fn, "tmp.gds",
    info.import=character(0), fmt.import=character(0))

# display
(f <- seqOpen("tmp.gds"))
seqClose(f)


# delete the temporary file
unlink("tmp.gds")
}

\keyword{gds}
\keyword{VCF}
\keyword{sequencing}
\keyword{genetics}
