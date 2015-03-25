SeqArray: Big Data Management of Genome-Wide Sequencing Variants
===


[![Build Status](https://travis-ci.org/zhengxwen/SeqArray.png)](https://travis-ci.org/zhengxwen/SeqArray)


## Features

Big data management of genome-wide variants using the CoreArray C++ library:
genotypic data and annotations are stored in an array-oriented manner,
offering efficient access of genetic variants using the S programming language.


## Bioconductor:

Development Version: v1.7.4

[http://www.bioconductor.org/packages/devel/bioc/html/SeqArray.html](http://www.bioconductor.org/packages/devel/bioc/html/SeqArray.html)


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
[gdsfmt](https://github.com/zhengxwen/gdsfmt/tarball/master),
[SeqArray](https://github.com/zhengxwen/SeqArray/tarball/master)
```sh
wget --no-check-certificate https://github.com/zhengxwen/gdsfmt/tarball/master -O gdsfmt_latest.tar.gz
wget --no-check-certificate https://github.com/zhengxwen/SeqArray/tarball/master -O SeqArray_latest.tar.gz
## Or
curl -L https://github.com/zhengxwen/gdsfmt/tarball/master/ -o gdsfmt_latest.tar.gz
curl -L https://github.com/zhengxwen/SeqArray/tarball/master/ -o SeqArray_latest.tar.gz

## Install
R CMD INSTALL gdsfmt_latest.tar.gz
R CMD INSTALL SeqArray_latest.tar.gz
```
