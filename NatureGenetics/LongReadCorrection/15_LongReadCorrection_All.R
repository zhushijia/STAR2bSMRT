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
folder="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/test/PackageReleaseResult5/"
system( paste("mkdir -p",folder) )


########################################################################################################################
##########################################   p3Del1: 581  ##############################################################
########################################################################################################################

LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24460_581/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R1.fastq"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R2.fastq"
outputDir=paste0(folder,"Case1_581")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )
    
########################################################################################################################
##########################################   p3Del2: 641  ##############################################################
########################################################################################################################

LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24461_641/polished_high_qv_consensus_isoforms.fasta"
#SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/641/641.R1.fastq"
#SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/641/641.R2.fastq"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/EF12_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/EF12_R2_001.fastq.gz"
outputDir=paste0(folder,"Case2_641")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )


########################################################################################################################
##########################################   Cont1: 2607  ##############################################################
########################################################################################################################

LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24463_2607/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R1.fastq"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R2.fastq"
outputDir=paste0(folder,"Cont1_2607")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )


########################################################################################################################
##########################################   Cont2: 553  ###############################################################
########################################################################################################################

LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24459_553/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/553/553.R1.fastq"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/553/553.R2.fastq"
outputDir=paste0(folder,"Cont2_553")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )



########################################################################################################################
##########################################   Cont3: 642  ###############################################################
########################################################################################################################


LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24462_642/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/642-2-redo_S14_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/642-2-redo_S14_R2_001.fastq.gz"
outputDir=paste0(folder,"Cont3_642")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )


########################################################################################################################
##########################################   Adult1: NRXN_adult_dlPFC1_12  #############################################
########################################################################################################################

LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26964/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/106781-adult-dIPFC_S3_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/106781-adult-dIPFC_S3_R2_001.fastq.gz"
outputDir=paste0(folder,"Adult1_NRXN_adult_dlPFC1_12")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )


########################################################################################################################
##########################################   Adult2: adult_dlPFC1_10  ##################################################
########################################################################################################################

LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26963/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/14-34-adult-dIPFC_S2_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/14-34-adult-dIPFC_S2_R2_001.fastq.gz"
outputDir=paste0(folder,"Adult2_adult_dlPFC1_10")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )


########################################################################################################################
##########################################   Adult3: adult_dlPFC1_13  ##################################################
########################################################################################################################

LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26965/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/12-27-adult-dIPFC_S7_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/12-27-adult-dIPFC_S7_R2_001.fastq.gz"
outputDir=paste0(folder,"Adult3_adult_dlPFC1_13")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )


########################################################################################################################
##########################################   Fetal1: fetal  ############################################################
########################################################################################################################

LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26962/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/fetal-diPFC_S1_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/fetal-diPFC_S1_R2_001.fastq.gz"
outputDir=paste0(folder,"Fetal1_fetal")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )


########################################################################################################################
##########################################   Fetal2: fetal 23 weeks  ###################################################
########################################################################################################################

LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027966/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/6_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/6_R2_001.fastq.gz"
outputDir=paste0(folder,"Fetal2_fetal_23wks")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )


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


########################################################################################################################
##########################################   GABA_553  #################################################################
########################################################################################################################

# 553 GABA
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027949/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/1_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/1_R2_001.fastq.gz"
outputDir=paste0(folder,"GABA_553")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )


########################################################################################################################
##########################################   NGN2_553  #################################################################
########################################################################################################################

# 553 NGN2
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027950/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/2_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/2_R2_001.fastq.gz"
outputDir=paste0(folder,"NGN2_553")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )


########################################################################################################################
##########################################   NGN2_641  #################################################################
########################################################################################################################

# 641 NGN2 
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027959/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/12_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/12_R2_001.fastq.gz"
outputDir=paste0(folder,"NGN2_641")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores )




# 581 GABA low number of phqv
#LR="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027964/polished_high_qv_consensus_isoforms.fasta"


# 581 NGN2
#LR="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027963/polished_high_qv_consensus_isoforms.fasta"


# 641 GABA
#LR="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/028024/polished_high_qv_consensus_isoforms.fasta"




