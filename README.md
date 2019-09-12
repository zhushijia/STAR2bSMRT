# STAR2bSMRT
STARlong and STARshort based Single Molecule Real Time Iso-Seq analysis

## Description
STAR2bSMRT is a novel hybrid sequencing-based alternative splicing identification method, which is specially tailored for long read RNAseq correction. 

The long read RNAseq differs from the long read DNAseq from three folds: 
1) as opposed to such primary application of DNAseq as the do novo genome assebmly, the RANseq-based splicing isoform identification often has the pre-built genome reference, which serves as the great prior knowledge for both read alignment and correction; 
2) among all sequences of the RNAseq reads, the very small profortion of splicing junction sites act the dominant role in splicing isoform identification; 
3) the RNAseq reads demonstrate variational coverages across the genome in response to different gene expression and alternative splicing. 

These difference motivates the novel long read RNAseq-specific correction method, STAR2bSMRT. Given the matched hybrid-sequencing, STAR2bSMRT first aligns both long and short reads to the genome, obtaining the approximate splicing junctions; next, differing from corretion of the whole long read suquence, STAR2bSMRT only corrects those junction sites, via maximizing the correlation between the long and short read junctions. 

STAR2bSMRT has the following contributions:
1) in addition to short reads, it also incorporates prior knowledge to perform long read correction: genome reference, and annotated junction sites (gtf file);
2) it only corrects junction sites, largely saving efforts from correcting the whole sequence;
3) it converts the question of sequence correction into that of statistical optimization, also enabling us to automatically select the best parameters. Two parameters are crucial during the optimization: thresSR (the number of short reads which support the splicing junction sites) and thresDis (the tolerance distance between short read and long read-derived junction sites). 

## Dependencies on packages
-  STARlong & STARshort - [STAR](https://github.com/alexdobin/STAR)
-  samtools - [samtools](http://samtools.sourceforge.net/)
-  bedtools - [bedtools](http://bedtools.readthedocs.io/en/latest/)
-  R package: [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html); [foreach](https://cran.r-project.org/web/packages/foreach/); [doMC](https://cran.r-project.org/web/packages/doMC/)


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

### 3) Install the STAR2bSMRT R package
```
git clone https://zhushijia:password@github.com/zhushijia/STAR2bSMRT.git 
R CMD INSTALL -l userFolder STAR2bSMRT
```
Alternatively, use [devtools](https://github.com/hadley/devtools) package
```
> install.packages("devtools")
> library(devtools)
> install_github("zhushijia/STAR2bSMRT")
```

## Tutorial
   See our [wiki](https://github.com/zhushijia/STAR2bSMRT/wiki)


