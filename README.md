# STAR2bSMRT
STARlong and STARshort based Single Molecule Real Time Iso-Seq analysis

## Description
We proposed a novel method for genome reference, annotation and short read based long read correction. A novel criterion was used to optimize the selection of thresholds for number of short reads supporting junction sites and difference between long read and short read detected junction sites. Since the correction only focus on the junction sites, avoiding from the large time concumption on the sequence correction, enabling ultrafast correction. 

## Dependencies on R packages
-  STARlong
-  STARshort
-  samtools
-  bedtools
-  Biostrings - [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
-  foreach - [foreach](https://cran.r-project.org/web/packages/foreach/)
-  doMC - [doMC](https://cran.r-project.org/web/packages/doMC/)


## Installation
### 1) Obtaining a recent version of R
Clear instructions for different version can be found here:
http://cran.fhcrc.org/

### 2) Install the dependent R packages
```
# install R package of Biostrings. 
# try http:// if https:// URLs are not supported
> source("https://bioconductor.org/biocLite.R")
> biocLite("Biostrings")

# install packages for parallel computating
> install.packages(c("foreach","doMC"))

```

### 3) Install the SMRTER R package
```
wget 
R CMD 
```
Alternatively, use [devtools](https://github.com/hadley/devtools) package
```
> install.packages("devtools")
> library(devtools)
> install_github("zhushijia/STAR2bSMRT")
```

## Tutorial
   See our [wiki](https://github.com/zhushijia/STAR2bSMRT/wiki)
   
  
