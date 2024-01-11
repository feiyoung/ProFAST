# ProFAST

=========================================================================
<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version-ago/ProFAST)](https://cran.r-project.org/package=ProFAST)
[![](https://cranlogs.r-pkg.org/badges/ProFAST?color=orange)](https://cran.r-project.org/package=ProFAST)
[![](https://cranlogs.r-pkg.org/badges/grand-total/ProFAST?color=orange)](https://cran.r-project.org/package=ProFAST)
<!-- badges: end -->

FAST is an advanced probabilistic factor analysis technique designed for spatially-aware dimension reduction in multi-section spatial transcriptomics data. The 'ProFAST' package incorporates the FAST method and is specifically developed by the Jin Liu's lab for the comprehensive analysis of multiple spatially resolved transcriptomics (SRT) datasets.

Check out  our [Package Website](https://feiyoung.github.io/ProFAST/index.html) for a more complete description of the methods and analyses. 

ProFAST  can be used to compare and contrast experimental datasets in a variety of contexts, for instance:

* Across experimental batches
* Across individuals
* Across different conditions (i.e., case and control)
* Across datasets with only partially shared cell/domain clusters

Once the embeddings of multiple datasets are estimated by ProFAST, the package provides functionality for further data exploration, 
analysis, and visualization. Users can:

* Align the embeddings obtained by ProFAST
* Identify clusters using all data information
* Recover comparable gene expression matrices among datasets
* Find significant shared (and dataset-specific) gene markers
* Visuzlize extracted embeddings using 3-dim tSNE and UMAP
* Visualize clusters and gene expression using tSNE and UMAP

# Installation
"ProFAST" depends on the 'Rcpp' and 'RcppArmadillo' package, which requires appropriate setup of computer. For the users that have set up system properly for compiling C++ files, the following installation command will work.
```{Rmd}
# Method 1: install ProFAST from CRAN
install.packages('ProFAST')


# For the newest version of ProFAST, users can use method 2 for installation.
# Method 2: Install ProFAST from Github
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("feiyoung/ProFAST")

# If some dependent packages (such as `scater`) on Bioconductor can not be installed nomrally, use following commands, then run abouve command.
if (!require("BiocManager", quietly = TRUE)) ## install BiocManager
    install.packages("BiocManager")
# install the package on Bioconducter
BiocManager::install(c("scater"))
```



## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Integration analysis for a toy example with three batches](https://feiyoung.github.io/ProFAST/articles/FASTsimu.html)
* [Integration  analysis for DLPFC data with two batches](https://feiyoung.github.io/ProFAST/articles/FASTdlpfc2.html)



For the users that don't have set up system properly, the following setup on different systems can be referred.
## Setup on Windows system
First, download [Rtools](https://cran.r-project.org/bin/windows/Rtools/); second, add the Rtools directory to the environment variable.


## Setup on MacOS system
First, install Xcode. Installation about Xcode can be referred [here](https://stackoverflow.com/questions/8291146/xcode-installation-on-mac).


Second, install "gfortran" for compiling C++ and Fortran at [here](https://github.com/fxcoudert/gfortran-for-macOS).


## Setup on Linux  system
If you use conda environment on Linux system and some dependent packages (such as `scater`) can not normally installed, you can search R package at anaconda.org website. We take the `scater` package as example, and its search result is https://anaconda.org/bioconda/bioconductor-scater. Then you can install it in conda environment by following command.
```{Linux}

conda install -c bioconda bioconductor-scater
```
For the user not using conda environment, if  dependent packages (such as `scater`) not normally installed are in Bioconductor, then use the following command to install the dependent packages.
```{Linux}
# install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# install the package on Bioconducter
BiocManager::install(c("scater"))
```
If  dependent packages (such as `DR.SC`) not normally installed are in CRAN, then use the following command to install the dependent packages.
```{Linux}
# install the package on CRAN
install.packages("DR.SC")
```


# Demonstration

For an example of typical ProFAST usage, please see our [Package Website](https://feiyoung.github.io/ProFAST/index.html) for a demonstration and overview of the functions included in ProFAST.

# NEWs
* ProFAST version 1.3 (2024-01-09)

* ProFAST version 1.1 (2023-06-03)


