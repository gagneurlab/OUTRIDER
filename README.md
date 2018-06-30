# OUTRIDER #
OUTRIDER is a tool to find aberrantly expressed genes in RNA-seq samples.
[Paper on bioRxiv](https://www.biorxiv.org/content/early/2018/06/14/322149)

[![Pipeline status](https://travis-ci.org/gagneurlab/OUTRIDER.svg?branch=master)](https://travis-ci.org/gagneurlab/OUTRIDER)
[![Version](https://img.shields.io/badge/Version-0.99.16-orange.svg)](https://github.com/gagneurlab/OUTRIDER/blob/master)
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


### Alternative tensorflow autoencoder

Additionally to the in R implemented autoencoder we provide a second 
experimental implementation written in Python using tensorflow.


#### Python version installation

For autoCorrection we need a running [Python3 version](https://www.python.org/downloads/)
and for OUTRIDER a running [R 3.4.x version](https://cran.r-project.org/).

It is a good practice to have a virtual environment in python to have a clean 
working environment and not to create conflicts with other projects. This is also
helpful for macOS users, due to a [matplotlib bug](https://matplotlib.org/faq/osx_framework.html#osxframework-faq).

```
python3 -m venv ~/python-env-OUTRIDER
source ~/python-env-OUTRIDER/bin/activate
```

#### autoCorrection and Python dependencies

The first part is to install the `autoCorrection` package written in python. We will use `pip` to
install the package and all dependencies. autoCorrection depends on TensorFlow, which can be installed
manualy if GPU acceleration or any other customization is wanted. For more details please see the 
TensorFlow instruction [here](https://www.tensorflow.org/install/). If the TensorFlow installation
is skipped by the user, TensorFlow will be installed automatically from pip.

Installation of autoCorrection package:

```
pip install --upgrade autoCorrection
```

#### Python Troubleshooting

##### Failed to create virtual env with anaconda

If you use anaconda as your python distribution `venv` does not work. Please use this workaround:

```
python3 -m venv ~/python-env-OUTRIDER --without-pip
source ~/python-env-OUTRIDER/bin/activate
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py
```

For more details, please have a look [here](https://stackoverflow.com/questions/38524856/anaconda-3-for-linux-has-no-ensurepip?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa)

##### Failed to import Tkinter

Tkinter should be installed with most python installations. If it is still missing, please install the extra package

```
sudo apt-get install python3-tk
```

