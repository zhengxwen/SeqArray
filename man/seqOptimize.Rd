\name{seqOptimize}
\alias{seqOptimize}
\title{Optimize the Storage of Data Array}
\description{
    Transpose data array or matrix for possibly higher-speed access.
}
\usage{
seqOptimize(gdsfn, target=c("chromosome", "by.sample"), format.var=TRUE,
    cleanup=TRUE, verbose=TRUE)
}
\arguments{
    \item{gdsfn}{the file name of GDS}
    \item{target}{"chromosome", "by.sample"; see details}
    \item{format.var}{a character vector for selected variable names,
        or \code{TRUE} for all variables, according to "annotation/format"}
    \item{cleanup}{call \code{link{cleanup.gds}} if \code{TRUE}}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    None.
}
\details{
    \code{"chromosome"}: adding or updating two additional nodes
    '@chrom_rle_val' and '@chrom_rle_len' for faster chromosome indexing,
    requiring SeqArray>=v1.20.0.

    \code{"by.sample"}: optimizing GDS file for
    \code{seqApply(..., margin="by.sample")}. Warning: optimizing GDS file for
    reading data by sample may increase file size by up to 2X as genotype data
    and all format data are duplicated.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqGetData}}, \code{\link{seqApply}}
}

\examples{
# the file name of VCF
(vcf.fn <- seqExampleFileName("vcf"))
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# convert
seqVCF2GDS(vcf.fn, "tmp.gds", storage.option="ZIP_RA")

# prepare data for the SeqVarTools package
seqOptimize("tmp.gds", target="by.sample")


# list the structure of GDS variables
(f <- seqOpen("tmp.gds"))
# close
seqClose(f)


# delete the temporary file
unlink("tmp.gds")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
