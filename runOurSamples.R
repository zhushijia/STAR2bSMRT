library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

genomeFasta = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
chrom = "chr2"
s = 50147488
e = 51259537
cores = 30
thresSR=c(1:100) 
thresDis=c(1:30)


########################################################################################################################
##########################################   Miseq  ####################################################################
########################################################################################################################

genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24460_581/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R2.fastq"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_adjustNCjunc/581"

STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc=TRUE , chrom , s , e , cores)
  



genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24461_641/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/641/641.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/641/641.R2.fastq"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/641"


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24463_2607/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R2.fastq"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/2607"


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24459_553/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/553/553.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/553/553.R2.fastq"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/553"


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24462_642/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/642-2-redo_S14_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/642-2-redo_S14_R2_001.fastq.gz"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/642"



########################################################################################################################
##########################################   Hiseq  ######################################################################
########################################################################################################################

genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24460_581/polished_high_qv_consensus_isoforms.fasta"
SR1="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_581-2-1-FBN/Sample_581-2-1-FBN.R1.fastq"
SR2="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_581-2-1-FBN/Sample_581-2-1-FBN.R2.fastq"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/HiSeq/581-2-1"


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24461_641/polished_high_qv_consensus_isoforms.fasta"
SR1="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_641-6-2-FBN/Sample_641-6-2-FBN.R1.fastq"
SR2="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_641-6-2-FBN/Sample_641-6-2-FBN.R2.fastq"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/HiSeq/641-6-2"


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24463_2607/polished_high_qv_consensus_isoforms.fasta"
SR1="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_2607-2-1-FBN/Sample_2607-2-1-FBN.R1.fastq"
SR2="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_2607-2-1-FBN/Sample_2607-2-1-FBN.R2.fastq"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/HiSeq/2607-2-1"


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24459_553/polished_high_qv_consensus_isoforms.fasta"
SR1="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_553-1-1-FBN/Sample_553-1-1-FBN.R1.fastq"
SR2="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_553-1-1-FBN/Sample_553-1-1-FBN.R1.fastq"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/HiSeq/553-1-1"



######################
genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26962/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/fetal-diPFC_S1_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/fetal-diPFC_S1_R2_001.fastq.gz"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/HumanBrain/fetal"


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26963/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/14-34-adult-dIPFC_S2_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/14-34-adult-dIPFC_S2_R2_001.fastq.gz"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/HumanBrain/adult_dlPFC1_10"


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26964/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/106781-adult-dIPFC_S3_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/106781-adult-dIPFC_S3_R2_001.fastq.gz"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/HumanBrain/NRXN_adult_dlPFC1_12"


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26965/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/12-27-adult-dIPFC_S7_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/12-27-adult-dIPFC_S7_R2_001.fastq.gz"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/HumanBrain/adult_dlPFC1_13"




