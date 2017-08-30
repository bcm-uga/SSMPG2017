# SSMPG 2017
Repository for the [Data challenge about Software and Statistical Methods for Population Genetics (SSMPG 2017)](https://data-institute.univ-grenoble-alpes.fr/data-institute/news-and-events/data-challenge-on-software-and-statistical-methods-for-population-genetics-ssmpg-2017--713800.htm) (Aussois, September 11-15 2017)

##  1. Install software

### Install R and R-studio
To participate to the challenge, you need to install [R](https://cran.r-project.org/) on your computer. To make R easier to use, we suggest to install [RStudio](https://www.rstudio.com/), which is an integrated development environment (IDE) for R.

### Install R packages (LEA/LFMM, OutFLANK, pcadapt, rehh)
To install R packages that are useful for the data challenge, copy and paste in R the following pieces of code

```r
#Install R packages for SSMPG 2017

#Package to install packages from github
require(devtools)

#Package to run LFMM
devtools::install_github("bcm-uga/LEA")

#Package to run OutFLANK
devtools::install_github("whitlock/OutFLANK")

#Package to run pcadapt
devtools::install_github("bcm-uga/pcadapt")

#Package q-value for controlling FDR
#Try https:// or http:// 
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")

#Package to run rehh
install.packages("rehh")

```

### Install BAYPASS, hapFLK, SelEstim, SweeD

**TO COMPLETE** by B Servin, R Vitalis, and P Pavlidis.

##  2. Download datasets

## 3. Form a team

## 4. Submit candidate markers 

## 5. Evaluation

The ranking of the teams will be based on the [G score](https://en.wikipedia.org/wiki/F1_score#G-measure). The G score varies from 0 (minimum score) to 1 (maximum score).

The G score depends on the false discovery rate (FDR), which is the percentage of false positive markers (or regions) in the submitted list, and of the power, which is the percentage of markers (or regions) involved in adaptation, which are found in the submitted list. Two G scores will be computed for each submission; one marker-based score evaluates whether submitted markers are correct and one region-based score evaluates whether submitted markers falls in correct 100 bp regions.

The mathematical definition of the G score is

![equation](https://latex.codecogs.com/gif.latex?G%20%3D%20%5Csqrt%7B%5Cmathrm%7B%281-FDR%29%7D%20%5Ccdot%20%5Cmathrm%7Bpower%7D%7D.)

For the 1st dataset, $G$ scores will be publicly available (public leaderboard). Participants should use the 1st dataset for training. **The final ranking of the teams will be based on the 2nd dataset for which $G$ scores are not public.** Prizes will be provided to the best team according to the marker-based score and to the best team according to the region-best score.
