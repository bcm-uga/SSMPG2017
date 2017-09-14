**New on github** 

Have a look to the [presentation](Presentations/Cichlid) that introduces the Cichlid challenge.

The final ranking (evaluation set) for the Dahu challenge is available [online](http://176.31.253.205/shiny/SSMPG2017/). **Congrats to all participants**

**Updated program**

Wednesday: Dahu challenge

* 12:00-16:00
Moutain walks

* 19:00
End of the Dahu challenge

* 19:00-19:30
Introduction about the Cichlid challenge

* 19:30-20:30
Dinner

* 20:30-22:00 Presentations about the Dahu challenge. Each presentation should contain 3 slides containing 
    1. Materials and Methods
    2. Results
    3. Discussion

Thursday: Cichlid challenge

* 10:00-10:30
Coffee break

* 12:00-13:30
Lunch

* 17:30-18:00
Coffee break

* 18:00-19:30 Presentations about the Cichlid challenge. Each presentation should contain 3 slides containing 
    1. Materials and Methods
    2. Results
    3. Discussion


* 20:00-21:30
Gala Dinner with [pierrade](https://i.pinimg.com/736x/3c/f3/46/3cf346cb13e8d0df5190b406388a7f68--fondue-et-raclette-sauces.jpg)


# SSMPG 2017
Repository for the [Data challenge about Software and Statistical Methods for Population Genetics (SSMPG 2017)](https://data-institute.univ-grenoble-alpes.fr/data-institute/news-and-events/data-challenge-on-software-and-statistical-methods-for-population-genetics-ssmpg-2017--713800.htm) (Aussois, September 11-15 2017)

[Introductory lecture about Genome Scans](https://docs.google.com/document/d/1wQjVDZ2ZTPZoGVjSvLdqv-IMfpuunxYKHysAop27Mlo/edit#)

##  1. Install software

### Install R and Rstudio
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
SweeD is hosted by the github of Nikolaos Alachiotis (https://github.com/alachins/sweed). It's a C software, so you need to download the source code and compile it in your machine. You need to:

- go to https://github.com/alachins/sweed

- Click on the Green button "Clone or download". You will see a link "https://github.com/alachins/sweed.git". 

- Copy that link (https://github.com/alachins/sweed.git)

- Go to your terminal and your preferred folder (e.g. software), and type:
``` 
git clone https://github.com/alachins/sweed.git
```

This will download the sweed source code from the github repository, and it will create a folder called sweed

- change directory within sweed

``` 
cd sweed
```
- There are several versions of the code. The simplest one uses only a single thread. To compile it you should type:
```
make -f Makefile.gcc
```

- If you want to compile the pthreads version, and be able to exploit multiple treads of your PC type:
```
make clean
make -f ./Makefile.PTHREADS.gcc
```
**YOu need to remove the \*.o files before compiling. The reason is that it's possible that \*.o files are associated to another version and you will not be able to produce eventually the executable**

Now, you should have a running version of SweeD within the sweed directory. 
To be able to run SweeD in any folder of your PC, please (i) add the sweed directory in your path, or (ii) copy the SweeD executable in a folder within your path. **Don't do both, there is no need**
For (i): open the .bashrc file with your favorite editor e.g. emacs, and add the line:
```
PATH=$PATH:<SWEED FOLDER>
export PATH
``` 
**Don't forget the PATH variable, or your PATH will not be correct and you will not be able to run many commands in Linux**
For (ii): assuming that $HOME/bin/ is in your PATH, simply:
``` 
cp SweeD $HOME/bin/
```

To run SweeD and assuming that the input file is called **input.vcf**: 

just type:
``` 
SweeD -input input.VCF -grid 1000 -name MYRUN
```
### Install OmegaPlus (NOT FOR THE CONTEST)
OmegaPlus detects the LD patterns *around* the target of beneficial mutation. To install:

```
git clone https://github.com/alachins/omegaplus.git
cd omegaplus
make clean
make
```

### Install RAiSD (NOT FOR THE CONTEST)
In contrast to SweeD, RAiSD uses all 3 signatures of a selective sweep, namely the reduction of polymorphism levels, the shift of the SFS, and the special patterns of LD. Instead of the full information for SFS and LD, it uses approximations based on SNP vectors. 

RAiSD is hosted by github at https://github.com/alachins/raisd

The following commands can be used to download and compile the source code. 

    $ mkdir RAiSD
    $ cd RAiSD
    $ wget https://github.com/alachins/raisd/archive/master.zip
    $ unzip master.zip
    $ cd raisd-master
    $ make
    
The executable is placed in the path RAiSD/raisd-master/bin/release. A link to the executable is placed in the installation folder, i.e., raisd-master.

##  2. Download datasets

The first challenge is the **Dahu** challenge. Teams are asked to analyze simulated data for the Dahu challenge.

Download the [training](https://drive.google.com/open?id=0B2RlEpeOlBn_VjF3ZnVocVpyeXM) and the [validation](https://drive.google.com/open?id=0B2RlEpeOlBn_NVlheHFTM2tmQjQ) dataset for the Dahu challenge.

The second challenge is the **Cichlid** challenge. Teams are asked to analyze true data for the  Cichlid challenge.

Download the [data](https://drive.google.com/open?id=0B2RlEpeOlBn_VXh3Rl9BSng5czA) for the Cichlid challenge.

## 3. Form a team

To participate to the challenge, you should form teams. A team can be composed of 1, 2, or 3 participants. Once you have chosen a team, click on "create teams" on the [submission website](http://176.31.253.205/shiny/SSMPG2017/).

## 4. Submit candidate markers 


The objective of the two data challenges is to find markers that are involved in adaptation.

To submit a list of markers involved in adaptation, you should use the [submission website](http://176.31.253.205/shiny/SSMPG2017/). The submitted file should be a .txt file with one column containing the indices of the candidate SNPs. An example of submission file containing a list of markers involved in adaptation is contained in the file [mysubmission.txt](https://drive.google.com/uc?export=download&id=0B9o4VIJJSodfOGhYUkFlalJ4NXM). 

Here is a piece of R code that shows how to make a submission file

```{r}
training<-readRDS("Presentations/pcadapt/sim1a.rds")
stat<-apply(training$G,FUN=mean,MARGIN=2)
write(ou<-order(stat,decreasing = FALSE)[1:100],"ridicule_never_killed_anyone.txt",ncolumns=1)
```

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

## 6. Simulation results
A comparison of results for the first QTL region. Vertical lines indicate quatitative trait nucleotides that have an effect on the trait under selection.
https://docs.google.com/document/d/1HTLXB-I4qWI6yp0u-ysqXq4224EWfJiVdT5bGayI63s/edit
