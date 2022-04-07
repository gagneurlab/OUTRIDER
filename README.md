# OUTRIDER2 #
OUTRIDER2 is a generalized framework for context-dependent outlier detection in -omics data.

The original method OUTRIDER is published in the [AJHG](https://doi.org/10.1016/j.ajhg.2018.10.025)
and available through [Bioconductor](http://bioconductor.org/packages/release/bioc/html/OUTRIDER.html).

[![Pipeline status](https://travis-ci.org/gagneurlab/OUTRIDER.svg?branch=outrider2)](https://travis-ci.org/gagneurlab/OUTRIDER)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/a2f6io5isq0jhobf/branch/outrider2?svg=true)](https://ci.appveyor.com/project/c-mertes/outrider/branch/outrider2)
[![Version](https://img.shields.io/badge/Version-1.99.0-green.svg)](https://github.com/gagneurlab/OUTRIDER/tree/outrider2)
[![Coverage status](https://codecov.io/gh/gagneurlab/OUTRIDER/branch/outrider2/graph/badge.svg)](https://codecov.io/github/gagneurlab/OUTRIDER?branch=outrider2)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/gagneurlab/OUTRIDER/blob/outrider2/LICENSE)


## OUTRIDER2 development installation

`OUTRIDER`(v2) is an R software package requiring a running [R 3.6 version or higher](https://cran.r-project.org/).

`OUTRIDER`(v2) is not yet on Bioconductor. Instead, it can be installed using `BiocManager` (incluidng its dependencies). 
This will use `remotes::install_github` to install the `OUTRIDER2` version from the `outrider2` branch on github. 

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("gagneurlab/OUTRIDER", ref="outrider2")
```

If you have dependency issues while installing any package, please have a look
at the Troubleshooting section or submit an issue on [GitHub](https://github.com/gagneurlab/OUTRIDER/issues).

### Alternative OUTRIDER2 development installation with devtools

If you want to get the latest development version of `OUTRIDER`(v2), you can also 
install it from GitHub with `devtools`. For this, you need a working development environment to compile the
C++ code (see for details: [Windows](https://cran.r-project.org/bin/windows/Rtools/)
or [MacOS X](https://cran.r-project.org/bin/macosx/tools/)).

```
install.packages('devtools')

# latest development version of OUTRIDER2 (on branch outrider2)
devtools::install_github('gagneurlab/OUTRIDER', ref='outrider2', dependencies=TRUE)

```

To check which versions/tags are available you can check the GitHub repo
[here](https://github.com/gagneurlab/OUTRIDER/releases).

### Quick tour through OUTRIDER

To get started with `OUTRIDER`, please have a look at our
[vignette](https://github.com/gagneurlab/OUTRIDER/tree/outrider2/vignettes/OUTRIDER2.pdf).
<!--In order to get the pdf version, please type the following code in an R session:

```
library(OUTRIDER)
vignette('OUTRIDER')
```
-->

### Toubleshooting

#### Missing libraries while compiling R packages

On some Linux distributions we need the developer libraries for compiling the R packages.

To install those packages, please run as administrator: 

For Ubuntu or debian based systems:
```
sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev zlib1g-dev libmysqld-dev
```

For centOS or RHEL based systems:
```
sudo yum install R-devel zlib-devel openssl-devel libcurl-devel libxml2-devel mariadb-devel
```
