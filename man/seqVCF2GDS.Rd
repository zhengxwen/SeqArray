\name{seqVCF2GDS}
\alias{seqVCF2GDS}
\alias{seqBCF2GDS}
\title{Reformat VCF Files}
\description{
    Reformats Variant Call Format (VCF) files.
}
\usage{
seqVCF2GDS(vcf.fn, out.fn, header=NULL, storage.option="LZMA_RA",
    info.import=NULL, fmt.import=NULL, genotype.var.name="GT",
    ignore.chr.prefix="chr", scenario=c("general", "imputation"),
    reference=NULL, start=1L, count=-1L, variant_count=NA_integer_,
    optimize=TRUE, raise.error=TRUE, digest=TRUE, use_Rsamtools=NA,
    parallel=FALSE, verbose=TRUE)
seqBCF2GDS(bcf.fn, out.fn, header=NULL, storage.option="LZMA_RA",
    info.import=NULL, fmt.import=NULL, genotype.var.name="GT",
    ignore.chr.prefix="chr", scenario=c("general", "imputation"),
    reference=NULL, optimize=TRUE, raise.error=TRUE, digest=TRUE,
    bcftools="bcftools", verbose=TRUE)
}
\arguments{
    \item{vcf.fn}{the file name(s) of VCF format; or a \code{\link{connection}}
        object}
    \item{bcf.fn}{a file name of binary VCF format (BCF)}
    \item{out.fn}{the file name of output GDS file}
    \item{header}{if NULL, \code{header} is set to be
        \code{\link{seqVCF_Header}(vcf.fn)}}
    \item{storage.option}{specify the storage and compression option,
        "ZIP_RA" (\code{\link{seqStorageOption}("ZIP_RA")}); or
        "LZMA_RA" to use LZMA compression algorithm with higher compression
        ratio by default; or "LZ4_RA" to use an extremely fast compression and
        decompression algorithm. "ZIP_RA.max", "LZMA_RA.max" and "LZ4_RA.max"
        correspond to the algorithms with a maximum compression level; the
        suffix "_RA" indicates that fine-level random access is available; see
        more details at \code{\link{seqStorageOption}}}
    \item{info.import}{characters, the variable name(s) in the INFO field
        for import; or \code{NULL} for all variables}
    \item{fmt.import}{characters, the variable name(s) in the FORMAT field
        for import; or \code{NULL} for all variables}
    \item{genotype.var.name}{the ID for genotypic data in the FORMAT column;
        "GT" by default (in VCF v4)}
    \item{ignore.chr.prefix}{a vector of character, indicating the prefix of
        chromosome which should be ignored, e.g., \code{"chr"}; it is not
        case-sensitive}
    \item{scenario}{"general": use float32 to store floating-point numbers
        (by default); "imputation": use packedreal16 to store DS and GP in the
        FORMAT field with four decimal place accuracy}
    \item{reference}{genome reference, like "hg19", "GRCh37"; if the genome
        reference is not available in VCF files, users could specify the
        reference here}
    \item{start}{the starting variant if importing part of VCF files}
    \item{count}{the maximum count of variant if importing part of VCF files,
        -1 indicates importing to the end}
    \item{variant_count}{\code{NA_integer_} (default) or a numeric vector
        specifying the numbers of variants in the VCF file(s) in \code{vcf.fn};
        only applicable when multiple cores are used; if the number of variants
        is known, the conversion can skip counting the variants before
        splitting the file(s); \code{variant_count} could be an approximate}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{raise.error}{\code{TRUE}: throw an error if numeric conversion fails;
        \code{FALSE}: get missing value if numeric conversion fails}
    \item{digest}{a logical value (TRUE/FALSE) or a character ("md5", "sha1",
        "sha256", "sha384" or "sha512"); add md5 hash codes to the GDS file
        if TRUE or a digest algorithm is specified}
    \item{use_Rsamtools}{only applicable when multiple cores are used;
        \code{use_Rsamtools} is passed to \code{\link{seqVCF_Header}} to get
        the total number of variants;
        \code{NA}: using Rsamtools when it is installed; \code{FALSE}: not use
        the Rsamtools package; \code{TRUE}: to use Rsamtools, if it is not
        installed, the function fails}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (parallel
        processing), a numeric value indicating the number of cores, or a
        cluster object for parallel processing; \code{parallel} is passed to
        the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details}
    \item{verbose}{if \code{TRUE}, show information}
    \item{bcftools}{the path of the program \code{bcftools}}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\details{
    If there are more than one files in \code{vcf.fn}, \code{seqVCF2GDS} will
merge all VCF files together if they contain the same samples. It is useful
to merge multiple VCF files if variant data are split by chromosomes.

    The real numbers in the VCF file(s) are stored in 32-bit floating-point
format by default. Users can set
\code{storage.option=seqStorageOption(float.mode="float64")}
to switch to 64-bit floating point format. Or packed real numbers can be
adopted by setting
\code{storage.option=seqStorageOption(float.mode="packedreal16:scale=0.0001")}.

    By default, the compression method is "LZMA_RA" (\url{https://tukaani.org/xz/},
LZMA algorithm with default compression level + independent data blocks for
fine-level random access). Users can maximize the compression ratio by
\code{storage.option="LZMA_RA.max"} or
\code{storage.option=seqStorageOption("LZMA_RA.max")}. LZMA is known to have
higher compression ratio than the zlib algorithm. LZ4
(\url{https://github.com/lz4/lz4}) is an option via
\code{storage.option="LZ4_RA"} or
\code{storage.option=seqStorageOption("LZ4_RA")}.

    If multiple cores/processes are specified in \code{parallel}, all VCF files
are scanned to calculate the total number of variants before format conversion,
and then split by the number of cores/processes.

    \code{storage.option="Ultra"} and \code{storage.option="UltraMax"} need much
larger memory than other compression methods. Users may consider using
\code{\link{seqRecompress}} to recompress the GDS file after calling
\code{seqVCF2GDS()} with \code{storage.option="ZIP_RA"}, since
\code{seqRecompress()} compresses data nodes one by one, taking much less memory
than \code{"Ultra"} and \code{"UltraMax"}.

    If \code{storage.option="LZMA_RA"} runs out of memory (e.g., there are too
many annotation fields in the VCF file), users could use
\code{storage.option="ZIP_RA"} and then call
\code{seqRecompress(, compress="LZMA")}.
}
\references{
    Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo,
M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., et al. (2011).
The variant call format and VCFtools. Bioinformatics 27, 2156-2158.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF_Header}}, \code{\link{seqStorageOption}},
    \code{\link{seqMerge}}, \code{\link{seqGDS2VCF}},
    \code{\link{seqRecompress}}
}

\examples{
# the VCF file
vcf.fn <- seqExampleFileName("vcf")

# conversion
seqVCF2GDS(vcf.fn, "tmp.gds", storage.option="ZIP_RA")

# conversion in parallel
seqVCF2GDS(vcf.fn, "tmp_p2.gds", storage.option="ZIP_RA", parallel=2L)


# display
(f <- seqOpen("tmp.gds"))
seqClose(f)



# convert without the INFO fields
seqVCF2GDS(vcf.fn, "tmp.gds", storage.option="ZIP_RA",
    info.import=character(0))

# display
(f <- seqOpen("tmp.gds"))
seqClose(f)



# convert without the INFO and FORMAT fields
seqVCF2GDS(vcf.fn, "tmp.gds", storage.option="ZIP_RA",
    info.import=character(0), fmt.import=character(0))

# display
(f <- seqOpen("tmp.gds"))
seqClose(f)


# delete the temporary file
unlink(c("tmp.gds", "tmp_p2.gds"), force=TRUE)
}

\keyword{gds}
\keyword{VCF}
\keyword{sequencing}
\keyword{genetics}
