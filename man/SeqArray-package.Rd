\name{SeqArray-package}
\alias{SeqArray-package}
\alias{SeqArray}
\docType{package}
\title{Big Data Management of Whole-Genome Sequence Variant Calls}
\description{
    Big-data management of whole-genome sequencing variants.
}
\details{
    As the cost of DNA sequencing rapidly decreases, whole-genome sequencing
(WGS) is generating data at an unprecedented rate. Scientists are being
challenged to manage data sets that are terabyte-sized, contain diverse types
of data and complex data relationships. Data analyses of WGS requires a general
file format for storing genetic variants including single nucleotide variations
(SNVs), insertions and deletions (indels) and structural variants. The variant
call format (VCF) is a generic and flexible format for storing DNA polymorphisms
developed for the 1000 Genomes Project that is the standard WGS format in use
today. VCF is a textual format usually stored in compressed files that supports
rich annotations and relatively efficient data retrieval. However, VCF files
are large and the computational burden associated with all data retrieval from
text files can be significant for a large WGS study with thousands of samples.

    To provide an efficient alternative to VCF for WGS data, we developed a new
data format and accompanying Bioconductor package, ``SeqArray''. Key features
of SeqArray are efficient storage including multiple high compression options,
data retrieval by variant or sample subsets, support for parallel access and
computing, and C++ integration in the R programming environment. The SeqArray
package provides R functions for efficient block-wise computations, and enables
scientists to develop custom R scripts for exploratory data analysis.

    Webpage: \url{http://github.com/zhengxwen/SeqArray},
        \url{http://www.bioconductor.org/packages/SeqArray/}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}

\examples{
# the file of VCF
vcf.fn <- seqExampleFileName("vcf")
vcf.fn
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# parse the header
seqVCF_Header(vcf.fn)

# get sample id
seqVCF_SampID(vcf.fn)

# convert
seqVCF2GDS(vcf.fn, "tmp.gds", storage.option="ZIP_RA")
seqSummary("tmp.gds")

# list the structure of GDS variables
f <- seqOpen("tmp.gds")
f

seqClose(f)
unlink("tmp.gds")


############################################################

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


# read multiple variables variant by variant
seqApply(f, c(geno="genotype", phase="phase", qual="annotation/id"),
    FUN=function(x) print(x), as.is="none")

# get the numbers of alleles per variant
head(seqApply(f, "allele",
    FUN=function(x) length(unlist(strsplit(x,","))), as.is="integer"))
# or
head(seqGetData(f, "$num_allele"))



################################################################

# remove the sample and variant filters
seqResetFilter(f)

# calculate the frequency of reference allele,
#   a faster version could be obtained by C coding
af <- seqApply(f, "genotype", FUN=function(x) mean(x==0L, na.rm=TRUE),
    as.is="double")
length(af)
summary(af)


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
