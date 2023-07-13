# [An R package `PICH`](https://github.com/hfang-bristol/PICH)

## @ Overview

> The `PICH` is an R package designed to identify therapeutic candidates (both therapeutic targets and repurposed drugs) based on multi-layered datasets related to clonal hematopoiesis (defined as the clonal expansion of the blood stem cell and its progeny).


## @ Installation

### 1. Install R

Please install R (version 4.2.1 or above); see https://cran.r-project.org

Below describes the instruction how to install R on `Ubuntu` (assuming you have a `ROOT (sudo)` privilege).

```ruby
sudo su
# here enter your password

wget http://www.stats.bris.ac.uk/R/src/base/R-4/R-4.2.1.tar.gz
tar xvfz R-4.2.1.tar.gz
cd ~/R-4.2.1
./configure
make
make check
make install
R # start R
```

### 2. Install R packages

```ruby
R # start R

# if the package 'BiocManager' not installed, please do so
if(!("BiocManager" %in% rownames(installed.packages()))) install.packages("BiocManager")

# first, install basic packages: remotes, tidyverse
BiocManager::install(c('remotes','tidyverse'), dependencies=T)

# then, install the package 'PICH' (now hosted at github)
BiocManager::install("hfang-bristol/PICH", dependencies=T, force=T)

# check the package 'PICH' successfully installed
library(PICH)
```


## @ Contact

> Please drop [email](mailto:fh12355@rjh.com.cn) for bug reports or enquiries.


