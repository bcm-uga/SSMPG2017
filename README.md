# SSMPG 2017
Repository for the [Data challenge about Software and Statistical Methods for Population Genetics (SSMPG 2017)](https://data-institute.univ-grenoble-alpes.fr/data-institute/news-and-events/data-challenge-on-software-and-statistical-methods-for-population-genetics-ssmpg-2017--713800.htm) (Aussois, September 11-15 2017)

[Introductory lecture about Genome Scans](https://docs.google.com/document/d/1wQjVDZ2ZTPZoGVjSvLdqv-IMfpuunxYKHysAop27Mlo/edit#)

##  1. Install software

### Install R and R-studio
To participate to the challenge, you need to install [R](https://cran.r-project.org/) on your computer. To make R easier to use, we suggest to install [RStudio](https://www.rstudio.com/), which is an integrated development environment (IDE) for R.

### Install R packages (LEA/LFMM, OutFLANK, pcadapt, rehh)
To install R packages that are useful for the data challenge, copy and paste in R the following pieces of code

```r
#Install R packages for SSMPG 2017

#Package to install packages from github
install.packages("devtools")

#Packages to run LEA and LFMM
devtools::install_github("bcm-uga/LEA")

install.packages("RSpectra")
devtools::install_github("bcm-uga/lfmm")

# Additional package 'cate' 
install.packages("cate")

#Package to run OutFLANK
devtools::install_github("whitlock/OutFLANK")

#Package to run pcadapt
devtools::install_github("bcm-uga/pcadapt")
install.packages("bigstatsr")
devtools::install_github("privefl/bigsnpr")

#Package q-value for controlling FDR
#Try https:// or http:// 
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")

#Package to run rehh
install.packages("rehh")

#Package to plot population trees
install.packages("ape")
```

### OutFLANK resources
[Link to tutorial for OutFLANK](http://rpubs.com/lotterhos/outflank)

[OutFLANK website](https://github.com/whitlock/OutFLANK)

[Whitlock MC and KE Lotterhos. 2015. Reliable Detection of Loci Responsible for Local Adaptation: Inference of a Null Model through Trimming the Distribution of FST. American Naturalist](http://www.journals.uchicago.edu/doi/abs/10.1086/682949)

### Install BAYPASS

Download the archive from http://www1.montpellier.inra.fr/CBGP/software/baypass/ or directly via the following command run on a terminal:
```
wget http://www1.montpellier.inra.fr/CBGP/software/baypass/files/baypass_2.1.tar.gz
```
Extract the archive, *e.g.*, from a terminal:
```
tar -zxvf baypass_2.1.tar.gz
```
The source files are to be found in the *src* subdirectory. BayPass is coded in Fortran90 and can therefore be compiled for any system supporting a Fortran90 compiler using the provided Makefile. This Makefile is designed to work with either the free compiler *gfortran* (if not already installed in your system, binaries are available at https://gcc.gnu.org/wiki/GFortranBinaries and are easy to install for most Windows, Mac and Linux OS versions) or the commercial *ifort* intel Fortran compiler. 
BayPass also uses OpenMP (http://openmp.org/wp/) to implement multithreading, which allows parallel calculation on computer systems that have multiple CPUs or CPUs with multiple cores. Users thus have to make sure that the corresponding libraries are installed (which is usually the case, on Linux OS or following compiler installation previously described). The following instructions run within the *src* subdirectory allows to compile the code and to produce a binary:
* using the *gfortran* free compiler (the command should automatically produce an executable called *g_baypass*):
```
make clean all FC=gfortran
```
* using the *ifort* intel Fortran compiler (the command should automatically produce an executable called *i_baypass*):
```
make clean all FC=ifort 
```
> Note: Under Linux (or MacOS), before the first use, make sure to give appropriate execution rights to the program. For instance you may run:
>```chmod +x baypass```

### Install hapFLK

hapflk is available as a python package. It has been tested to work on Linux and MacOSX. Before installing hapflk, you will need to install [python 2.7](https://www.python.org/downloads/) and [numpy and scipy](https://www.scipy.org/install.html). You also need a C compiler (e.g. gcc) but this should be the case already. Once this is done, hapflk can be installed using pip (copy paste the following in a terminal):

```
sudo pip install hapflk
```

In the future, hapflk can be upgraded using :

```
sudo pip install hapflk --upgrade
```

Checkout the [hapflk webpage](https://forge-dga.jouy.inra.fr/projects/hapflk/)
for some documentation and companion scripts.

### Install SelEstim

Download the archive from http://www1.montpellier.inra.fr/CBGP/software/selestim/, or using the following command line from a terminal:
```
wget http://www1.montpellier.inra.fr/CBGP/software/selestim/files/SelEstim_1.1.7.zip
```
Extract the archive, e.g., from a terminal:
```
unzip SelEstim_1.1.7.zip
```
The source files are to be found in the src/ subdirectory of that archive. SelEstim is coded using C programming language and can therefore be compiled for any system supported by [gcc](http://gcc.gnu.org/). To do so, Windows users may need to get a [gcc](http://gcc.gnu.org/), e.g. by installing [MinGW](http://www.mingw.org/), [mingw-64](http://mingw-w64.org/doku.php), or [Cygwin](https://sourceware.org/cygwin/). To compile the code and get the selestim binary, use the provided Makefile in the src/ subdirectory:
```
make clean all
```
> Note: with Linux (or Mac OS), before the first use, make sure to give appropriate execution rights to the program. For instance you may run:
>```chmod +x selestim```

SelEstim uses [OpenMP](href{http://openmp.org/wp/) to implement multithreading, which allows parallel calculation on on computer systems that have multiple CPUs or CPUs with multiple cores. Make sure that the corresponding libraries are installed, which is typically the case on Linux, Mac OS and Windows (provided the above recommendations for installation of gcc have been followed). 
> Note: The [gcc](http://gcc.gnu.org/) version included with OS X may generate executable code that results in runtime error (Abort trap: 6) when more than one thread is used. In that case, you first need to install a recent version of [gcc](http://gcc.gnu.org/), following the instructions at http://hpc.sourceforge.net/. Then, you can recompile SelEstim using the following instruction:
> ```make clean all CC=/usr/local/bin/gcc```
> (assuming gcc has been installed in the /usr/local/ subdirectory.)

### Install SweeD

**TO COMPLETE** by P Pavlidis.

##  2. Download datasets

The first challenge is the **Dahu** challenge. Teams are asked to analyze simulated data for the Dahu challenge.

Download the [training](https://drive.google.com/open?id=0B2RlEpeOlBn_Rk52bEgwZmtzeGc) and the [validation](https://drive.google.com/open?id=0B2RlEpeOlBn_RndaUExoQTFkaUU) dataset for the Dahu challenge.

The second challenge is the **Cichlid** challenge. Teams are asked to analyze true data for the  Cichlid challenge.

Download the [data](https://drive.google.com/open?id=0B2RlEpeOlBn_SC03ZVBVY1VCMW8) for the Cichlid challenge.

## 3. Form a team

To participate to the challenge, you should form teams. A team can be composed of 1, 2, or 3 participants. Once you have chosen a team, click on "create teams" on the [submission website](http://176.31.253.205/shiny/SSMPG2017/).

## 4. Submit candidate markers 


The objective of the two data challenges is to find markers that are involved in adaptation.

To submit a list of markers involved in adaptation, you should use the [submission website](http://176.31.253.205/shiny/SSMPG2017/). The submitted file should be a .txt file with one column containing the indices of the candidate SNPs. An example of submission file containing a list of markers involved in adaptation is contained in the file [mysubmission.txt](https://drive.google.com/uc?export=download&id=0B9o4VIJJSodfOGhYUkFlalJ4NXM). 

## 5. Evaluation

### Dahu challenge
Submissions will be evaluated by comparing submitted list to the list of causal adaptive SNPs.

The ranking of the teams will be based on the [G score](https://en.wikipedia.org/wiki/F1_score#G-measure). The G score varies from 0 (minimum score) to 1 (maximum score).

The G score depends on the false discovery rate (FDR), which is the percentage of false positive markers (or regions) in the submitted list, and of the power, which is the percentage of markers (or regions) involved in adaptation, which are found in the submitted list. Two G scores will be computed for each submission; one marker-based score evaluates whether submitted markers are correct and one region-based score evaluates whether submitted markers fall in correct 100 bp regions.

The mathematical definition of the G score is

![](https://latex.codecogs.com/gif.latex?G%20%3D%20%5Csqrt%7B%5Cmathrm%7B%281-FDR%29%7D%20%5Ccdot%20%5Cmathrm%7Bpower%7D%7D.)

For the *training* dataset, G scores will be publicly available (public leaderboard). Participants should use the *training* dataset for training. 

**The final ranking of the teams will be based on the  *validation*  dataset for which G scores are not public.** 

A prize will be provided to the best team according to the marker-based score and another prize will be provided to the best team according to the region-best score.

### Cichlid challenge

Submissions will be evaluated based on subjective evaluations by instructors based on presentations.

## 5. Presentations

During the SSMPG prize ceremony, each team will be asked to present 2-3 slides for each challenge.

