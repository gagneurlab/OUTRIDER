# OUTRIDER #
OUTRIDER is a tool to find aberrantly expressed genes in RNA-seq samples.
The method is published in the [AJHG](https://doi.org/10.1016/j.ajhg.2018.10.025)
and available through [Bioconductor](http://bioconductor.org/packages/release/bioc/html/OUTRIDER.html).

[![Pipeline status](https://travis-ci.org/gagneurlab/OUTRIDER.svg?branch=master)](https://travis-ci.org/gagneurlab/OUTRIDER)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/a2f6io5isq0jhobf/branch/master?svg=true)](https://ci.appveyor.com/project/c-mertes/outrider/branch/master)
[![Version](https://img.shields.io/badge/Version-1.5.1-green.svg)](https://github.com/gagneurlab/OUTRIDER/tree/master)
[![Coverage status](https://codecov.io/gh/gagneurlab/OUTRIDER/branch/master/graph/badge.svg)](https://codecov.io/github/gagneurlab/OUTRIDER?branch=master)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/gagneurlab/OUTRIDER/blob/master/LICENSE)


## Installation

`OUTRIDER` is an R software package requiring a running [R 3.6 version or higher](https://cran.r-project.org/).

We will use `BiocManager` to install the package and its dependencies. If you
want to install the devel version of `OUTRIDER`, please install it 
with devtools from GitHub directly [(see below)](#OUTRIDER-development-installation).


```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install('OUTRIDER')
```

If you have dependency issues while installing any package, please have a look
at the Troubleshooting section or submit an issue on [GitHub](https://github.com/gagneurlab/OUTRIDER/issues).

### OUTRIDER development installation

If you want to get the latest
development version of `OUTRIDER`, you can install it from GitHub with
`devtools`. For this, you need a working development environment to compile the
C++ code (see for details: [Windows](https://cran.r-project.org/bin/windows/Rtools/)
or [MacOS X](https://cran.r-project.org/bin/macosx/tools/)).

```
install.packages('devtools')

# latest development version
devtools::install_github('gagneurlab/OUTRIDER', dependencies=TRUE)

# installing a specific version/tag of OUTRIDER
devtools::install_github('gagneurlab/OUTRIDER@1.3.5', dependencies=TRUE)
```

To check which versions/tags are available you can check the GitHub repo
[here](https://github.com/gagneurlab/OUTRIDER/releases).

### Quick tour through OUTRIDER

To get started with `OUTRIDER`, please have a look at our
[vignette](http://bioconductor.org/packages/release/bioc/vignettes/OUTRIDER/inst/doc/OUTRIDER.pdf).
In order to get the pdf version, please type the following code in an R session:

```
library(OUTRIDER)
vignette('OUTRIDER')
```

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
