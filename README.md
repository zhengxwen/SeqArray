SeqArray: Big Data Management of Genome-wide Sequencing Variants
===


[![Build Status](https://travis-ci.org/zhengxwen/SeqArray.png)](https://travis-ci.org/zhengxwen/SeqArray)


## Features

Big data management of genome-wide variants using the CoreArray library: genotypic data and annotations are stored in an array-oriented manner, offering efficient access of genetic variants using the R language.


## Installation

* Bioconductor repository:
```
source("http://bioconductor.org/biocLite.R")
library(BiocInstaller)
BiocInstaller::useDevel()

biocLite("SeqArray")
```

* Development version from Github (v0.99.1):
```
library("devtools")
install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SeqArray")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually.
