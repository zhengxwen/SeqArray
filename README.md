SeqArray: Big Data Management of Whole-genome Sequence Variant Calls
===

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)

[![Availability](http://www.bioconductor.org/shields/availability/release/SeqArray.svg)](http://www.bioconductor.org/packages/release/bioc/html/SeqArray.html)
[![Years-in-BioC](http://www.bioconductor.org/shields/years-in-bioc/SeqArray.svg)](http://www.bioconductor.org/packages/release/bioc/html/SeqArray.html)
[![Build Status](https://travis-ci.org/zhengxwen/SeqArray.png)](https://travis-ci.org/zhengxwen/SeqArray)
[![Build status](https://ci.appveyor.com/api/projects/status/noil0942el3iohqs?svg=true)](https://ci.appveyor.com/project/zhengxwen/seqarray)
[![codecov.io](https://codecov.io/github/zhengxwen/SeqArray/coverage.svg?branch=master)](https://codecov.io/github/zhengxwen/SeqArray?branch=master)


## Features

Big data management of whole-genome sequence variant calls with thousands of individuals: genotypic data (e.g., SNVs, indels and structural variation calls) and annotations in SeqArray files are stored in an array-oriented and compressed manner, with efficient data access using the R programming language.

The SeqArray package is built on top of Genomic Data Structure (GDS) data format, and defines required data structure for a SeqArray file. GDS is a flexible and portable data container with hierarchical structure to store multiple scalable array-oriented data sets. It is suited for large-scale datasets, especially for data which are much larger than the available random-access memory. It also offers the efficient operations specifically designed for integers of less than 8 bits, since a diploid genotype usually occupies fewer bits than a byte. Data compression and decompression are available with relatively efficient random access. A high-level R interface to GDS files is available in the package gdsfmt (http://bioconductor.org/packages/gdsfmt).


## Bioconductor:

Release Version: v1.12.9

[http://www.bioconductor.org/packages/release/bioc/html/SeqArray.html](http://www.bioconductor.org/packages/release/bioc/html/SeqArray.html)

* [Help Documents](http://zhengxwen.github.io/SeqArray/release/help/00Index.html)
* Tutorials: [Data Management](http://www.bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/SeqArrayTutorial.html), [R Integration](http://www.bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/R_Integration.html), [Overview Slides](http://www.bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/OverviewSlides.html)

Development Version: v1.13.5

[http://www.bioconductor.org/packages/devel/bioc/html/SeqArray.html](http://www.bioconductor.org/packages/devel/bioc/html/SeqArray.html)

* [Help Documents](http://zhengxwen.github.io/SeqArray/devel/help/00Index.html)
* Tutorials: [Data Management](http://www.bioconductor.org/packages/devel/bioc/vignettes/SeqArray/inst/doc/SeqArrayTutorial.html), [R Integration](http://www.bioconductor.org/packages/devel/bioc/vignettes/SeqArray/inst/doc/R_Integration.html), [Overview Slides](http://www.bioconductor.org/packages/devel/bioc/vignettes/SeqArray/inst/doc/OverviewSlides.html)



## Installation (requiring >=R_v3.3.0)

* Bioconductor repository:
```R
source("http://bioconductor.org/biocLite.R")
biocLite("SeqArray")
```

* Development version from Github:
```R
library("devtools")
install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SeqArray")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](http://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.

* Install the package from the source code:
[gdsfmt](https://github.com/zhengxwen/gdsfmt), [SeqArray](https://github.com/zhengxwen/SeqArray)
```sh
wget --no-check-certificate https://github.com/zhengxwen/gdsfmt/tarball/master -O gdsfmt_latest.tar.gz
wget --no-check-certificate https://github.com/zhengxwen/SeqArray/tarball/master -O SeqArray_latest.tar.gz
R CMD INSTALL gdsfmt_latest.tar.gz
R CMD INSTALL SeqArray_latest.tar.gz

## Or
curl -L https://github.com/zhengxwen/gdsfmt/tarball/master/ -o gdsfmt_latest.tar.gz
curl -L https://github.com/zhengxwen/SeqArray/tarball/master/ -o SeqArray_latest.tar.gz
R CMD INSTALL gdsfmt_latest.tar.gz
R CMD INSTALL SeqArray_latest.tar.gz
```



## SeqArray File Download

* [1000 Genomes Project](http://bochet.gcc.biostat.washington.edu/seqarray/1000genomes)



## Examples

```R
library(SeqArray)

gds.fn <- seqExampleFileName("gds")

# open a GDS file
f <- seqOpen(gds.fn)

# display the contents of the GDS file
f

# close the file
seqClose(f)
```

```R
## Object of class "SeqVarGDSClass"
## File: SeqArray/extdata/CEU_Exon.gds (387.3K)
## +    [  ] *
## |--+ description   [  ] *
## |--+ sample.id   { VStr8 90 ZIP_ra(30.8%), 222B }
## |--+ variant.id   { Int32 1348 ZIP_ra(35.7%), 1.9K }
## |--+ position   { Int32 1348 ZIP_ra(86.4%), 4.6K }
## |--+ chromosome   { VStr8 1348 ZIP_ra(2.66%), 91B }
## |--+ allele   { VStr8 1348 ZIP_ra(17.2%), 928B }
## |--+ genotype   [  ] *
## |  |--+ data   { Bit2 2x90x1348 ZIP_ra(28.4%), 16.8K } *
## |  |--+ ~data   { Bit2 2x1348x90 ZIP_ra(36.0%), 21.3K } *
## |  |--+ extra.index   { Int32 3x0 ZIP_ra, 17B } *
## |  \--+ extra   { Int16 0 ZIP_ra, 17B }
## |--+ phase   [  ]
## |  |--+ data   { Bit1 90x1348 ZIP_ra(0.36%), 55B } *
## |  |--+ ~data   { Bit1 1348x90 ZIP_ra(0.36%), 55B } *
## |  |--+ extra.index   { Int32 3x0 ZIP_ra, 17B } *
## |  \--+ extra   { Bit1 0 ZIP_ra, 17B }
## |--+ annotation   [  ]
## |  |--+ id   { VStr8 1348 ZIP_ra(41.0%), 5.8K }
## |  |--+ qual   { Float32 1348 ZIP_ra(0.91%), 49B }
## |  |--+ filter   { Int32,factor 1348 ZIP_ra(0.89%), 48B } *
## |  |--+ info   [  ]
## |  |  |--+ AA   { VStr8 1348 ZIP_ra(24.2%), 653B } *
## |  |  |--+ AC   { Int32 1348 ZIP_ra(27.2%), 1.4K } *
## |  |  |--+ AN   { Int32 1348 ZIP_ra(20.6%), 1.1K } *
## |  |  |--+ DP   { Int32 1348 ZIP_ra(62.6%), 3.3K } *
## |  |  |--+ HM2   { Bit1 1348 ZIP_ra(117.2%), 198B } *
## |  |  |--+ HM3   { Bit1 1348 ZIP_ra(117.2%), 198B } *
## |  |  |--+ OR   { VStr8 1348 ZIP_ra(14.0%), 238B } *
## |  |  |--+ GP   { VStr8 1348 ZIP_ra(34.4%), 5.3K } *
## |  |  \--+ BN   { Int32 1348 ZIP_ra(21.6%), 1.1K } *
## |  \--+ format   [  ]
## |     \--+ DP   [  ] *
## |        |--+ data   { Int32 90x1348 ZIP_ra(33.8%), 160.3K }
## |        \--+ ~data   { Int32 1348x90 ZIP_ra(32.2%), 152.8K }
## \--+ sample.annotation   [  ]
##    \--+ family   { VStr8 90 ZIP_ra(34.7%), 135B }
```


## Key Functions in the SeqArray Package

| Function     | Description |
|:-------------|:-------------------------------------------|
| seqVCF2GDS   | Reformat VCF files. [![](vignettes/link.png)](http://zhengxwen.github.io/SeqArray/release/help/seqVCF2GDS.html)  |
| seqSetFilter | Define a data subset of samples or variants. [![](vignettes/link.png)](http://zhengxwen.github.io/SeqArray/release/help/seqSetFilter.html)  |
| seqGetData   | Get data from a SeqArray file with a defined filter. [![](vignettes/link.png)](http://zhengxwen.github.io/SeqArray/release/help/seqGetData.html)  |
| seqApply     | Apply a user-defined function over array margins. [![](vignettes/link.png)](http://zhengxwen.github.io/SeqArray/release/help/seqApply.html)  |
| seqParallel  | Apply functions in parallel. [![](vignettes/link.png)](http://zhengxwen.github.io/SeqArray/release/help/seqParallel.html)  |


## Significant User-visible Changes (since v1.11.16)

* `seqSummary(gds, "genotype")$seldim` returns a vector with 3 integers (ploidy, # of selected samples, # of selected variants) instead of 2 integers

