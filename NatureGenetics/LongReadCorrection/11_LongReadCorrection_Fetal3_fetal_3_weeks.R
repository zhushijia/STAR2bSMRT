#STARlong_PATH="/sc/orga/projects/schzrnas/sjzhu/bitbucket/STAR/STAR/source/"
#STARshort_PATH="/sc/orga/projects/schzrnas/sjzhu/bitbucket/STARshort/STAR/source/"
#samtools_PATH="/hpc/packages/minerva-common/samtools/1.1/bin/"
#bedtools_PATH="/hpc/packages/minerva-common/BEDTools/2.27.1/bin/"
#kallisto_PATH="/hpc/packages/minerva-common/kallisto/0.45.0/kallisto_linux-v0.45.0/"
#export PATH=$PATH:$STARlong_PATH:$STARshort_PATH:$samtools_PATH:$bedtools_PATH:$kallisto_PATH


library(STAR2bSMRT,lib.loc="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/forRelease/r2")

genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
genomeFasta = "/hpc/users/xzhus01/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
chrom = "chr2"
s = 50147488
e = 51259537
cores = 10
thresSR=c(1:100) 
thresDis=c(0:30)
adjustNCjunc=TRUE
fixedMatchedLS=FALSE
useSJout=FALSE
fuzzyMatch=100
folder="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/test/PackageReleaseResult2/"
system( paste("mkdir -p",folder) )



########################################################################################################################
##########################################   Fetal3: fetal 3 weeks  ####################################################
########################################################################################################################

LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027957/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/7_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/7_R2_001.fastq.gz"
outputDir=paste0(folder,"Fetal3_fetal_3wks")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )

