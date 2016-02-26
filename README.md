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

Big data management of whole-genome sequence variant calls with thousands of individuals: genotypic data (e.g., SNVs, indels and structural variation calls) and annotations in GDS files are stored in an array-oriented and compressed manner, with efficient data access using the R programming language.

## Bioconductor:

Release Version: v1.10.6

[http://www.bioconductor.org/packages/release/bioc/html/SeqArray.html](http://www.bioconductor.org/packages/release/bioc/html/SeqArray.html)

* [Help Documents](http://zhengxwen.github.io/SeqArray/release/help/00Index.html)
* Tutorials: [Data Management](http://www.bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/SeqArrayTutorial.html), [Data Analytics](http://www.bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/AnalysisTutorial.html)

Development Version: v1.11.13

[http://www.bioconductor.org/packages/devel/bioc/html/SeqArray.html](http://www.bioconductor.org/packages/devel/bioc/html/SeqArray.html)

* [Help Documents](http://zhengxwen.github.io/SeqArray/devel/help/00Index.html)
* Tutorials: [Data Management](http://www.bioconductor.org/packages/devel/bioc/vignettes/SeqArray/inst/doc/SeqArrayTutorial.html), [Data Analytics](http://www.bioconductor.org/packages/devel/bioc/vignettes/SeqArray/inst/doc/AnalysisTutorial.html), [Overview Slides](http://www.bioconductor.org/packages/devel/bioc/vignettes/SeqArray/inst/doc/OverviewSlides.html)



## Installation

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
## |--+ sample.id   { VStr8 90 ZIP_RA(30.83%), 222B }
## |--+ variant.id   { Int32 1348 ZIP_RA(35.72%), 1.9K }
## |--+ position   { Int32 1348 ZIP_RA(86.44%), 4.6K }
## |--+ chromosome   { VStr8 1348 ZIP_RA(2.66%), 91B }
## |--+ allele   { VStr8 1348 ZIP_RA(17.19%), 928B }
## |--+ genotype   [  ] *
## |  |--+ data   { Bit2 2x90x1348 ZIP_RA(28.39%), 16.8K } *
## |  |--+ ~data   { Bit2 2x1348x90 ZIP_RA(36.04%), 21.3K } *
## |  |--+ extra.index   { Int32 3x0 ZIP_RA, 17B } *
## |  |--+ extra   { Int16 0 ZIP_RA, 17B }
## |--+ phase   [  ]
## |  |--+ data   { Bit1 90x1348 ZIP_RA(0.36%), 55B } *
## |  |--+ ~data   { Bit1 1348x90 ZIP_RA(0.36%), 55B } *
## |  |--+ extra.index   { Int32 3x0 ZIP_RA, 17B } *
## |  |--+ extra   { Bit1 0 ZIP_RA, 17B }
## |--+ annotation   [  ]
## |  |--+ id   { VStr8 1348 ZIP_RA(41.02%), 5.8K }
## |  |--+ qual   { Float32 1348 ZIP_RA(0.91%), 49B }
## |  |--+ filter   { Int32,factor 1348 ZIP_RA(0.89%), 48B } *
## |  |--+ info   [  ]
## |  |  |--+ AA   { VStr8 1348 ZIP_RA(24.22%), 653B } *
## |  |  |--+ AC   { Int32 1348 ZIP_RA(27.23%), 1.4K } *
## |  |  |--+ AN   { Int32 1348 ZIP_RA(20.62%), 1.1K } *
## |  |  |--+ DP   { Int32 1348 ZIP_RA(62.57%), 3.3K } *
## |  |  |--+ HM2   { Bit1 1348 ZIP_RA(117.16%), 198B } *
## |  |  |--+ HM3   { Bit1 1348 ZIP_RA(117.16%), 198B } *
## |  |  |--+ OR   { VStr8 1348 ZIP_RA(13.98%), 238B } *
## |  |  |--+ GP   { VStr8 1348 ZIP_RA(34.36%), 5.3K } *
## |  |  |--+ BN   { Int32 1348 ZIP_RA(21.64%), 1.1K } *
## |  |--+ format   [  ]
## |  |  |--+ DP   [  ] *
## |  |  |  |--+ data   { Int32 90x1348 ZIP_RA(33.83%), 160.3K }
## |  |  |  |--+ ~data   { Int32 1348x90 ZIP_RA(32.23%), 152.8K }
## |--+ sample.annotation   [  ]
## |  |--+ family   { VStr8 90 ZIP_RA(34.70%), 135B }
```
