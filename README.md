# OUTRIDER #

[![Pipeline status](https://travis-ci.org/gagneurlab/OUTRIDER.svg?branch=master)](https://travis-ci.org/gagneurlab/OUTRIDER)
[![Coverage status](https://codecov.io/gh/gagneurlab/OUTRIDER/branch/master/graph/badge.svg)](https://codecov.io/github/gagneurlab/OUTRIDER?branch=master)
[![Coverage status](https://codecov.io/gh/c-mertes/OUTRIDER/branch/master/graph/badge.svg)](https://codecov.io/github/c-mertes/OUTRIDER?branch=master)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/gagneurlab/OUTRIDER/blob/master/LICENSE)

## Installation

`OUTRIDER` consists of two software packages. Therefore both requirements have 
to be met before you can run the full software and pipeline.

### Prerequisite

For autoCorrection we need a running [Python3 version](https://www.python.org/downloads/)
and for OUTRIDER a running [R 3.4.x version](https://cran.r-project.org/).

It is a good practice to have a virtual environment in python to have a clean 
working environment and not to start conflicts with other projects. This is also
helpful for macOS users, due to a [matplotlib bug](https://matplotlib.org/faq/osx_framework.html#osxframework-faq).

```
python -m venv ~/python-env-OUTRIDER
source ~/python-env-OUTRIDER/bin/activate
```

For OUTRIDER we will us `devtools` and `BiocInstaller` to install it and its dependencies.

```
Rscript -e "install.packages('devtools', repos='http://cran.us.r-project.org')"
Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('BiocInstaller')"
```

If you have dependency issues while installing devtools please have a look at the Troubleshooting section.

### autoCorrection and Python dependencies

The first part is to install the `autoCorrection` package written in python. We will use `pip` to
install the package and all dependencies. autoCorrection depends on TensorFlow, which can be installed
manualy if GPU acceleration or any other customization is wanted. For more details please see the 
TensorFlow instruction [here](https://www.tensorflow.org/install/). If the tensorFlow installation
is skipped by the user, tensorFlow will be installed automatically from pip.

Installation of autoCorrection package:

```
pip install --upgrade autoCorrection
```

### OUTRIDER and R dependencies

With autoCorrection installed on our system, we can install the OUTRIDER package in R.

```
Rscript -e 'devtools::install_github('gagneurlab/OUTRIDER', dependencies=TRUE, build_vignettes=TRUE)'
``` 

### Quick tour through OUTRIDER

To get startet with OUTRIDER please have a look at our [vignette](vignette/OUTRIDER.Rnw)
In order to get the pdf version please type the following code in an R sessions:

```
library(OUTRIDER)
vignette('OUTRIDER')
```

### Toubleshooting

#### Missing libraries while install R packages

On some Linux distributions we need the developer libraries for compiling the R packages.

To install those packages please run as adriministrator: 

For Ubuntu or debian based systems:
```
sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev zlib1g-dev libmysqld-dev
```

For centOS or RHEL based systems:
```
sudo yum install R-devel zlib-devel openssl-devel libcurl-devel libxml2-devel mariadb-devel
```

#### Failed to create virtual env with anaconda

If you use anaconda as your python distribution `venv` does not work. Please use this workaround:

```
python -m venv ~/python-env-OUTRIDER --without-pip
source ~/python-env-OUTRIDER/bin/activate
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py
```

For more details please have a look [here](https://stackoverflow.com/questions/38524856/anaconda-3-for-linux-has-no-ensurepip?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa)
