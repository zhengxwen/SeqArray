\name{SeqArray-package}
\alias{SeqArray-package}
\alias{SeqArray}
\docType{package}
\title{Big Data Management of Genome-wide Sequencing Variants}
\description{
    Big-data management of genome-wide sequencing variants.
}
\details{
    In the era of big data, thousands of gigabyte-size data sets are
challenging scientists for data management, even on well-equipped hardware.
Currently, next-generation sequencing techniques are being adopted to
investigate common and rare variants, making the analyses of large-scale
genotypic data challenging. For example, the 1000 Genomes Project has
identified approximately 38 million single nucleotide polymorphisms (SNPs),
1.4 million short insertions and deletions, and more than 14,000 larger
deletions from whole-genome sequencing technologies. In the near future, new
technologies, like third-generation whole-genome sequencing, will be enabling
data to be generated at an unprecedented scale. The Variant Call Format (VCF)
was developed for the 1000 Genomes Project, which is a generic text format for
storing DNA polymorphism data such as SNPs, insertions, deletions and
structural variants, together with rich annotations. However, this format is
less efficient for large-scale analyses since numeric data have to be parsed
from a text VCF file before further analyses. The computational burden
associated with sequencing variants is especially evident with large sample
and variant sizes, and it requires efficient numerical implementation and
data management.

    Here I introduce a high-performance C/C++ computing library CoreArray
(\url{http://corearray.sourceforge.net}) for big-data management of genome-wide
variants. CoreArray was designed for developing portable and scalable storage
technologies for bioinformatics data, allowing parallel computing at the
multicore and cluster levels. It provides the genomic data structure (GDS) file
format for array-oriented data: this is a universal data format to store
multiple data variables in a single file. A hierarchical data structure is used
to store multiple extensible data variables in the GDS format, and all datasets
are stored in a single file with chunked storage layout. Here, I focus on the
application of CoreArray for statisticians working in the R environment, and
developed an R/Bioconductor package SeqArray to address or reduce the
computational burden associated with data management of sequencing variants.
The kernels of SeqArray are written in C/C++ and highly optimized. Genotypic
data and annotations are stored in an array-oriented manner, offering efficient
access of genetic variants using the R language. There are five key functions
in SeqArray, and most of data analyses could be done using these 6 functions:
\tabular{ll}{
    \bold{Function} \tab \bold{Description} \cr
    seqVCF2GDS \tab Imports VCF files \cr
    seqSummary \tab Gets the summary of a sequencing GDS file (# of samples,
        # of variants, INFO/FORMAT variables, etc) \cr
    seqSetFilter \tab Sets a filter to sample or variant (define a subset
        of data) \cr
    seqGetData \tab Gets data from a sequencing GDS file (from a subset
        of data) \cr
    seqApply \tab Applies a user-defined function over array margins \cr
    seqParallel \tab Applies functions in parallel \cr
}

The 1000 Genomes Project released 39 million genetic variants for 1092
individuals, and a 26G data file was created by SeqArray to store sequencing
variants with phasing information, where 2 bits were used as an atomic data
type. The file size can be further reduced to 1.3G by compression algorithms
without sacrificing access efficiency, since it has a large proportion of
rare variants.

SeqArray will be of great interest to scientists involved in data analyses of
large-scale genomic sequencing data using R environment, particularly those
with limited experience of low-level C programming and parallel computing.

    Webpage: \url{http://corearray.sourceforge.net/},
        \url{http://github.com/zhengxwen/SeqArray}

    Tutorial: \url{http://corearray.sourceforge.net/tutorials/SeqArray/}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}

\examples{
# the file of VCF
vcf.fn <- seqExampleFileName("vcf")
vcf.fn
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# parse the header
seqVCF.Header(vcf.fn)

# get sample id
seqVCF.SampID(vcf.fn)

# convert
seqVCF2GDS(vcf.fn, "tmp.gds")
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
seqApply(f, "allele",
    FUN=function(x) length(unlist(strsplit(x,","))), as.is="integer")


################################################################

# remove the sample and variant filters
seqSetFilter(f)

# calculate the frequency of reference allele,
#   a faster version could be obtained by C coding
af <- seqApply(f, "genotype", FUN=function(x) mean(x==0, na.rm=TRUE), as.is="double")
length(af)
summary(af)


################################################################
# run in parallel

library(parallel)

# Use option cl.core to choose an appropriate cluster size or number of cores
cl <- makeCluster(getOption("cl.cores", 2))

# run in parallel
afreq <- seqParallel(cl, f, FUN = function(gdsfile) {
        seqApply(gdsfile, "genotype", as.is="double",
            FUN=function(x) mean(x==0, na.rm=TRUE))
    }, split = "by.variant")

length(afreq)
summary(afreq)


stopCluster(cl)


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
