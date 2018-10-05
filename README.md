# OUTRIDER #
OUTRIDER is a tool to find aberrantly expressed genes in RNA-seq samples.
[Paper on bioRxiv](https://www.biorxiv.org/content/early/2018/06/14/322149)

[![Pipeline status](https://travis-ci.org/gagneurlab/OUTRIDER.svg?branch=master)](https://travis-ci.org/gagneurlab/OUTRIDER)
[![Version](https://img.shields.io/badge/Version-0.99.28-orange.svg)](https://github.com/gagneurlab/OUTRIDER/tree/master)
[![Coverage status](https://codecov.io/gh/gagneurlab/OUTRIDER/branch/master/graph/badge.svg)](https://codecov.io/github/gagneurlab/OUTRIDER?branch=master)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/gagneurlab/OUTRIDER/blob/master/LICENSE)

## Installation

`OUTRIDER` is an R software package requiring a running [R 3.4 version or higher](https://cran.r-project.org/).

### Prerequisite

We will use `devtools` and `BiocInstaller` to install the package and its dependencies.

```
install.packages('devtools', repos='http://cran.us.r-project.org')
source('https://bioconductor.org/biocLite.R')
biocLite('BiocInstaller')
```

If you have dependency issues while installing devtools, please have a look at the Troubleshooting section.


### OUTRIDER and R dependencies

The `OUTRIDER` R package and its dependencies can then be installed by running:

```
devtools::install_github('gagneurlab/OUTRIDER', dependencies=TRUE)
``` 

### Quick tour through OUTRIDER

To get started with `OUTRIDER`, please have a look at our [vignette](https://i12g-gagneurweb.in.tum.de/public/paper/OUTRIDER/OUTRIDER-vignette.pdf).
In order to get the pdf version, please type the following code in an R session:

```
library(OUTRIDER)
vignette('OUTRIDER')
```

### Toubleshooting

#### Missing libraries while install R packages

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