\name{seqGDS2VCF}
\alias{seqGDS2VCF}
\title{Convert to a VCF file}
\description{
    Convert a GDS file to a VCF file.
}
\usage{
seqGDS2VCF(gdsfile, vcf.fn, info.var=NULL, fmt.var=NULL, verbose=TRUE)
}
\arguments{
    \item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
    \item{vcf.fn}{the file name, output a file of VCF format}
    \item{info.var}{a list of variable names in the INFO field, or NULL for
        using all variables; \code{character(0)} for no variable
        in the INFO field}
    \item{fmt.var}{a list of variable names in the FORMAT field, or NULL for
        using all variables; \code{character(0)} for no variable
        in the FORMAT field}
    \item{verbose}{if TRUE, show information}
}
\value{
    Return the file name of VCF file with an absolute path.
}
\details{
    \code{\link{seqSetFilter}} can be used to define a subset of data for
the export.

    GDS -- Genomic Data Structures used for storing genetic array-oriented
        data, and the file format used in the \link{gdsfmt} package.

    VCF -- The Variant Call Format (VCF), which is a generic format for
        storing DNA polymorphism data such as SNPs, insertions, deletions and
        structural variants, together with rich annotations.
}
\references{
    The variant call format and VCFtools.
    Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE,
    Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project Analysis Group.
    Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.

    \url{http://corearray.sourceforge.net/}
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqVCF2GDS}}
}

\examples{
# the GDS file
(gds.fn <- seqExampleFileName("gds"))

# display
(f <- seqOpen(gds.fn))

# output the first 10 samples
samp.id <- seqGetData(f, "sample.id")
seqSetFilter(f, sample.id=samp.id[1:5])


# convert
seqGDS2VCF(f, "tmp.vcf.gz")

# no INFO and FORMAT
seqGDS2VCF(f, "tmp1.vcf.gz", info.var=character(0), fmt.var=character(0))

# output BN,GP,AA,DP,HM2 in INFO (the variables are in this order), no FORMAT
seqGDS2VCF(f, "tmp2.vcf.gz", info.var=c("BN","GP","AA","DP","HM2"), fmt.var=character(0))


# read
(txt <- readLines("tmp.vcf.gz", n=20))
(txt <- readLines("tmp1.vcf.gz", n=20))
(txt <- readLines("tmp2.vcf.gz", n=20))





#########################################################################
# Users could compare the new VCF file with the original VCF file
# call "diff" in Unix (a command line tool comparing files line by line)

# using all samples and variants
seqSetFilter(f)

# convert
seqGDS2VCF(f, "tmp.vcf.gz")


# file.copy(seqExampleFileName("vcf"), "old.vcf.gz", overwrite=TRUE)
# system("diff <(gunzip -c old.vcf.gz) <(gunzip -c tmp.vcf.gz)")

# 1a2,3
# > ##fileDate=20130309
# > ##source=SeqArray_RPackage_v1.0

# LOOK GOOD!


# delete temporary files
unlink(c("tmp.vcf.gz", "tmp1.vcf.gz", "tmp2.vcf.gz"))

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{VCF}
\keyword{sequencing}
\keyword{genetics}
