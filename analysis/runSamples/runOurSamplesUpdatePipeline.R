library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone4/setup")

adjustNCjunc=FALSE
useSJout=FALSE
fuzzyMatch=0

for ( adjustNCjunc in c(FALSE,TRUE) )
{
  for ( useSJout in c(FALSE,TRUE) )
  {
    for ( fuzzyMatch in c(0,100) )
    {
      testParameters(adjustNCjunc , fixedMatchedLS , useSJout , fuzzyMatch )
    }
  }
}


testParameters = function(adjustNCjunc , fixedMatchedLS , useSJout , fuzzyMatch )
{
  
  genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
  genomeFasta = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
  chrom = "chr2"
  s = 50147488
  e = 51259537
  cores = 30
  thresSR=c(1:100) 
  thresDis=c(1:30)
  #adjustNCjunc=FALSE
  fixedMatchedLS=FALSE
  #useSJout=FALSE
  #fuzzyMatch=0
  folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/"
  system( paste("mkdir -p",folder) )
  
  
  ########################################################################################################################
  ##########################################   Miseq  ####################################################################
  ########################################################################################################################
  
  LRphqv="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24461_641/polished_high_qv_consensus_isoforms.fasta"
  #SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/641/641.R1.fastq"
  #SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/641/641.R2.fastq"
  SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/EF12_R1_001.fastq.gz"
  SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/EF12_R2_001.fastq.gz"
  outputDir=paste0(folder,"641")
  #SoutputDir = paste0(outputDir,"/SR")
  SoutputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping/EF12"
  LoutputDir = paste0(outputDir,"/LR")
  STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                   SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                   thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                   fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                   chrom=chrom , s=s , e=e , cores=cores ,
                   SoutputDir=SoutputDir , LoutputDir=LoutputDir )
  
  
  
  LRphqv="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24463_2607/polished_high_qv_consensus_isoforms.fasta"
  SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R1.fastq"
  SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R2.fastq"
  outputDir=paste0(folder,"2607")
  SoutputDir = paste0(outputDir,"/SR")
  LoutputDir = paste0(outputDir,"/LR")
  STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                   SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                   thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                   fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                   chrom=chrom , s=s , e=e , cores=cores ,
                   SoutputDir=SoutputDir , LoutputDir=LoutputDir )
  
  
  LRphqv="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24460_581/polished_high_qv_consensus_isoforms.fasta"
  SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R1.fastq"
  SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R2.fastq"
  outputDir=paste0(folder,"581")
  SoutputDir = paste0(outputDir,"/SR")
  LoutputDir = paste0(outputDir,"/LR")
  STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                   SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                   thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                   fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                   chrom=chrom , s=s , e=e , cores=cores ,
                   SoutputDir=SoutputDir , LoutputDir=LoutputDir )
  
  
  LRphqv="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24459_553/polished_high_qv_consensus_isoforms.fasta"
  SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/553/553.R1.fastq"
  SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/553/553.R2.fastq"
  outputDir=paste0(folder,"553")
  SoutputDir = paste0(outputDir,"/SR")
  LoutputDir = paste0(outputDir,"/LR")
  STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                   SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                   thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                   fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                   chrom=chrom , s=s , e=e , cores=cores ,
                   SoutputDir=SoutputDir , LoutputDir=LoutputDir )
  
  LRphqv="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24462_642/polished_high_qv_consensus_isoforms.fasta"
  SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/642-2-redo_S14_R1_001.fastq.gz"
  SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/642-2-redo_S14_R2_001.fastq.gz"
  outputDir=paste0(folder,"642")
  SoutputDir = paste0(outputDir,"/SR")
  LoutputDir = paste0(outputDir,"/LR")
  STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                   SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                   thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                   fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                   chrom=chrom , s=s , e=e , cores=cores ,
                   SoutputDir=SoutputDir , LoutputDir=LoutputDir )
  
  
  #LRphqv="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24461_641/polished_high_qv_consensus_isoforms.fasta"
  #SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/641-FB-neurons_S6_R1_001.fastq.gz"
  #SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/641-FB-neurons_S6_R2_001.fastq.gz"
  #outputDir=paste0(folder,"641_SR2")
  #STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)
  
  ########################################################################################################################
  ##########################################   hipsc and adult  ####################################################################
  ########################################################################################################################
  
  LRphqv="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26962/polished_high_qv_consensus_isoforms.fasta"
  SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/fetal-diPFC_S1_R1_001.fastq.gz"
  SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/fetal-diPFC_S1_R2_001.fastq.gz"
  outputDir=paste0(folder,"fetal")
  SoutputDir = paste0(outputDir,"/SR")
  LoutputDir = paste0(outputDir,"/LR")
  STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                   SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                   thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                   fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                   chrom=chrom , s=s , e=e , cores=cores ,
                   SoutputDir=SoutputDir , LoutputDir=LoutputDir )
  
  LRphqv="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26963/polished_high_qv_consensus_isoforms.fasta"
  SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/14-34-adult-dIPFC_S2_R1_001.fastq.gz"
  SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/14-34-adult-dIPFC_S2_R2_001.fastq.gz"
  outputDir=paste0(folder,"adult_dlPFC1_10")
  SoutputDir = paste0(outputDir,"/SR")
  LoutputDir = paste0(outputDir,"/LR")
  STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                   SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                   thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                   fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                   chrom=chrom , s=s , e=e , cores=cores ,
                   SoutputDir=SoutputDir , LoutputDir=LoutputDir )
  
  LRphqv="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26964/polished_high_qv_consensus_isoforms.fasta"
  SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/106781-adult-dIPFC_S3_R1_001.fastq.gz"
  SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/106781-adult-dIPFC_S3_R2_001.fastq.gz"
  outputDir=paste0(folder,"NRXN_adult_dlPFC1_12")
  SoutputDir = paste0(outputDir,"/SR")
  LoutputDir = paste0(outputDir,"/LR")
  STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                   SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                   thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                   fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                   chrom=chrom , s=s , e=e , cores=cores ,
                   SoutputDir=SoutputDir , LoutputDir=LoutputDir )
  
  LRphqv="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26965/polished_high_qv_consensus_isoforms.fasta"
  SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/12-27-adult-dIPFC_S7_R1_001.fastq.gz"
  SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/12-27-adult-dIPFC_S7_R2_001.fastq.gz"
  outputDir=paste0(folder,"adult_dlPFC1_13")
  SoutputDir = paste0(outputDir,"/SR")
  LoutputDir = paste0(outputDir,"/LR")
  STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                   SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                   thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                   fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                   chrom=chrom , s=s , e=e , cores=cores ,
                   SoutputDir=SoutputDir , LoutputDir=LoutputDir )

}


