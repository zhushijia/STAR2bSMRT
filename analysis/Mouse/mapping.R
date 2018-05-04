library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
library(Biostrings)

genomeDir="/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/gencode/GRCm38.p5/genome/StarIndexUsingPrimaryGtf"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/SRR1184043_4.fa"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/mouse-dIPFC-redo_S13_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/mouse-dIPFC-redo_S13_R2_001.fastq.gz"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT"

LoutputDir = paste0(outputDir,"/LR")
SoutputDir = paste0(outputDir,"/SR")
EoutputDir = paste0(outputDir,"/Exp")

system( paste0( "mkdir -p " , LoutputDir ) )
system( paste0( "mkdir -p " , SoutputDir ) )
system( paste0( "mkdir -p " , EoutputDir ) )

starLong( genomeDir , LR , LoutputDir )
starShort( genomeDir , SR1 , SR2 , SoutputDir )

