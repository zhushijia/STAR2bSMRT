folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"





export PATH=$PATH:/sc/orga/projects/schzrnas/sjzhu/bitbucket/STAR/STAR/source/
export PATH=$PATH:/sc/orga/projects/schzrnas/sjzhu/bitbucket/STARshort/STAR/source/
export PATH=$PATH:/hpc/packages/minerva-common/samtools/1.1/bin/
export PATH=$PATH:/hpc/packages/minerva-common/BEDTools/2.27.1/bin/

library(STAR2bSMRT,lib.loc="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/forRelease/r1")

genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
genomeFasta = "/hpc/users/xzhus01/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
chrom = "chr2"
s = 50147488
e = 51259537
cores = 30
thresSR=c(1:100) 
thresDis=c(0:30)
adjustNCjunc=TRUE
fixedMatchedLS=FALSE
useSJout=FALSE
fuzzyMatch=100
folder="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/test/PackageReleaseResult/"
system( paste("mkdir -p",folder) )


###########################################
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24461_641/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/EF12_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/EF12_R2_001.fastq.gz"
 
outputDir=paste0(folder,"641_2")


###########################################
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24463_2607/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R1.fastq"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R2.fastq"
outputDir=paste0(folder,"2607_2")
SoutputDir = paste0(outputDir,"/SR")
LoutputDir = paste0(outputDir,"/LR")
  

  