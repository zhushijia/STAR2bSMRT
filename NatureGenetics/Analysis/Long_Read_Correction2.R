export PATH=$PATH:/sc/orga/projects/schzrnas/sjzhu/bitbucket/STAR/STAR/source/
export PATH=$PATH:/sc/orga/projects/schzrnas/sjzhu/bitbucket/STARshort/STAR/source/
export PATH=$PATH:/hpc/packages/minerva-common/samtools/1.1/bin/
export PATH=$PATH:/hpc/packages/minerva-common/BEDTools/2.27.1/bin/


library(STAR2bSMRT,lib.loc="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

genomeFasta = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
chrom = "chr2"
s = 50147488
e = 51259537
cores = 30
thresSR=c(1:100) 
thresDis=c(1:30)
adjustNCjunc=FALSE
fixedMatchedLS=FALSE
folder="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/"
system( paste("mkdir -p",folder) )

########################################################################################################################
##########################################   Miseq  ####################################################################
########################################################################################################################

genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24461_641/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/641-FB-neurons_S6_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/641-FB-neurons_S6_R2_001.fastq.gz"
outputDir=paste0(folder,"641_SR2")
STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24460_581/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R2.fastq"
outputDir=paste0(folder,"581")
STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)

genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24461_641/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/641/641.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/641/641.R2.fastq"
outputDir=paste0(folder,"641")
STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)



genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24463_2607/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R2.fastq"
outputDir=paste0(folder,"2607")
STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24459_553/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/553/553.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/553/553.R2.fastq"
outputDir=paste0(folder,"553")
STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24462_642/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/642-2-redo_S14_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/642-2-redo_S14_R2_001.fastq.gz"
outputDir=paste0(folder,"642")
STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)


######################
genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26962/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/fetal-diPFC_S1_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/fetal-diPFC_S1_R2_001.fastq.gz"
outputDir=paste0(folder,"fetal")
STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26963/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/14-34-adult-dIPFC_S2_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/14-34-adult-dIPFC_S2_R2_001.fastq.gz"
outputDir=paste0(folder,"adult_dlPFC1_10")
STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26964/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/106781-adult-dIPFC_S3_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/106781-adult-dIPFC_S3_R2_001.fastq.gz"
outputDir=paste0(folder,"NRXN_adult_dlPFC1_12")
STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_primerR2/26965/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/12-27-adult-dIPFC_S7_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/12-27-adult-dIPFC_S7_R2_001.fastq.gz"
outputDir=paste0(folder,"adult_dlPFC1_13")
STAR2bSMRT( genomeDir , genomeFasta , LR , SR1 , SR2 , thresSR , thresDis , outputDir , adjustNCjunc , fixedMatchedLS , chrom , s , e , cores)


library(STAR2bSMRT,lib.loc="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone4/setup")

STAR2bSMRT_NRXN <- function( genomeDir, genomeFasta, LRphqv=NULL, LRflnc=NULL, LRnfl=NULL,
                             SR1, SR2=NULL, useSJout=TRUE,  adjustNCjunc=FALSE, 
                             thresSR, thresDis, outputDir, fixedMatchedLS=FALSE, fuzzyMatch=100, 
                             chrom=NULL , s=0 , e=Inf , cores=10 ,
                             SoutputDir , LoutputDir )
{
    
    library(Biostrings)
    library(foreach)
    library(doMC)
    registerDoMC(cores)
    
    
    ############################################################################
    ############   STARshort mapping and junction sites for short reads
    ############################################################################
    
    #SoutputDir = paste0(outputDir,"/SR")
    SRalignment = paste0(SoutputDir,"/alignments.bam")
    system( paste0( "mkdir -p " , SoutputDir ) )
    #starShort( genomeDir , SR1 , SR2 , SoutputDir )
    
    if( useSJout )
    {
        SRjunc = getJuncBySJout( SJout="SJ.out.tab", SoutputDir, chrom=chrom, s=s, e=e )
    } else {
        SRjunc = getJuncBySam( SRalignment, SoutputDir, chrom=chrom, s=s, e=e )
    }
    
    
    ############################################################################
    ############   STARlong mapping and junction sites for long reads
    ############################################################################
    
    ############   for phqv	############   
    
    if( !is.null(LRphqv)  )
    {
        #LoutputDir = paste0(outputDir,"/LR")
        LRalignment = paste0(LoutputDir,"/Aligned.out.sam")
        system( paste0( "mkdir -p " , LoutputDir ) )
        #starLong( genomeDir=genomeDir , LR=LRphqv , outputDir=LoutputDir , cores=cores , SJ=NULL )
        
        LRread = getReadByJI( LRalignment , LoutputDir )
        LRread = subset(LRread , start < 50150000 )
        exp = phqvExp(LRphqv,LoutputDir)  # get coverage for all phqv 
        LRread = merge( LRread , exp , by="id" ) 
        LRread$coverage = LRread$full_length_coverage
        LRinfo = getLRinfo( LRread ,  chrom=chrom , s=s , e=e )
        LRread = LRinfo$LRread
        LRjunc = LRinfo$LRjunc
        LRtag = LRinfo$LRtag
    }
    
    ############   for flnc	############
    
    if( is.null(LRphqv) & !is.null(LRflnc) & is.null(LRnfl) )
    {
        LoutputDir = paste0(outputDir, "/LR_flnc")
        LRalignment = paste0(LoutputDir, "/Aligned.out.sam")
        system( paste0( "mkdir -p " , LoutputDir ) )
        starLong( genomeDir=genomeDir, LR=LRflnc, outputDir=LoutputDir, cores=cores, SJ=NULL )
        
        LRread = getReadByJI( LRalignment, LoutputDir )
        LRread$group = sapply( strsplit(as.character(LRread$id),"/") , function(x) x[1] )
        LRinfo = getLRinfo( LRread,  chrom=chrom, s=s, e=e )
        LRread = LRinfo$LRread
        LRjunc = LRinfo$LRjunc
        LRtag = LRinfo$LRtag
    }
    
    ############   for both flnc and nfl	############
    
    if( is.null(LRphqv) & !is.null(LRflnc) & !is.null(LRnfl) )
    {
        LoutputDir1 = paste0(outputDir, "/LR_flnc")
        LRalignment1 = paste0(LoutputDir1, "/Aligned.out.sam")
        system( paste0( "mkdir -p " , LoutputDir1 ) )
        starLong( genomeDir=genomeDir, LR=LRflnc, outputDir=LoutputDir1, cores=cores, SJ=NULL )
        
        LoutputDir2 = paste0(outputDir, "/LR_nfl")
        LRalignment2 = paste0(LoutputDir2, "/Aligned.out.sam")
        system( paste0( "mkdir -p " , LoutputDir2 ) )
        starLong( genomeDir=genomeDir, LR=LRnfl, outputDir=LoutputDir2, cores=cores, SJ=NULL )
        
        
        LRread1 = getReadByJI( LRalignment1, LoutputDir1 )
        LRread1$group = sapply( strsplit(as.character(LRread1$id),"/") , function(x) x[1] )
        LRread2 = getReadByJI( LRalignment1, LoutputDir2 )
        LRread2$group = sapply( strsplit(as.character(LRread2$id),"/") , function(x) x[1] )
        LRread = rbind( LRread1 , LRread2 )
        
        LRinfo1 = getLRinfo( LRread1,  chrom=chrom, s=s, e=e )
        LRinfo = getLRinfo( LRread,  chrom=chrom, s=s, e=e )
        
        LRread = LRinfo1$LRread # LRread for flnc
        LRtag = LRinfo1$LRtag # LRtag for flnc 
        LRjunc = LRinfo$LRjunc # LRjunc for both flnc and nfl
    }
    
    ############################################################################
    ############   grid searching
    ############################################################################
    
    if( fixedMatchedLS )
    {
        matchedLS = matchLSjunc( LRjunc , SRjunc )
    } else {
        matchedLS = NULL
    }
    
    score = gridSearch( LRjunc , SRjunc , thresSR , thresDis , adjustNCjunc , matchedLS , fuzzyMatch )
    
    ij = which( score==max(score) , arr.ind=T )
    ts = thresSR[ ij[1,1] ]
    td = thresDis[ ij[1,2] ]
    cat( ts , td , score[ij] , '\n ')
    
    correction = generateCorrectedIsoform( LRjunc , SRjunc, LRtag , LRread  , ts , td , matchedLS , fuzzyMatch )
    print(correction[[1]][c(2,3)])
    
    
    EoutputDir = paste0(folder,basename(outputDir),"/adjustNCjunc_",adjustNCjunc,
                        "_fixedMatchedLS_",fixedMatchedLS,"_useSJout_",useSJout,
                        "_fuzzyMatch_",fuzzyMatch)
    
    system( paste0( "mkdir -p " , EoutputDir ) )
    
    setwd( EoutputDir )
    genome = readDNAStringSet(genomeFasta)
    
    ###############################################################################################################
    
    pdf( paste0( "JuncExp_LR_ts",ts,"_td",td,".pdf") )
    
    juncExp = do.call( rbind, lapply( correction , function(x) x$LSjuncCount ))
    lrCount = log10(juncExp$lrCount)
    srCount = log10(juncExp$srCount)
    juncCorr = cor.test(srCount,lrCount,method='spearman')$estimate
    cols = sapply( juncExp$motif , function(x) ifelse(x==0,1,2) )
    cols[juncExp$motif==1] = 3
    plot( lrCount , srCount , col=cols , pch=17 , main=paste0("JuncExp by Long and Short Reads: r=",signif(juncCorr,3)) ,  xlab="Log10 Long Read" , ylab="Log10 Short Read"  )
    abline(lm( srCount~lrCount ))
    
    par(mfrow=c(2,1))
    log10fc = lrCount - srCount
    JuncNames = paste(juncExp$start , juncExp$end)
    barplot( log10fc , cex.names=0.6 , col=cols , ylab="log10(lrCount/srCount)", names=JuncNames , las=3 )
    
    dev.off()
    
    ###############################################################################################################
    tag = paste0( "exp", correction[[chrom]]$normalizedIsoformCount )
    exonList = juncToExon( juncList=correction[[chrom]]$isoform , s=50149082 , e=51255411 , tag=tag )
    #writeGff( isoform=correction[[chrom]]$isoform , file = gffName , exp=correction[[chrom]]$exp , chrom='chr2' , s=50149082 , e=51255411 )
    
    ###############################################################################################################
    #seq = generateSeq( genome=genome , isoform=correction[[chrom]]$isoform , exp=correction[[chrom]]$exp , chrom='chr2' , s=50149082 , e=51255411  )
    seq = generateSeq( genome=genome , isoform=exonList )
    fastaName = paste0( "isoform_ts",ts,"_td",td,".fa")
    translated = sapply( seq$translated , function(x) ifelse( x , "translated" , "untranslated" )  )
    names(seq$dna) = paste(names(seq$dna),translated,sep="_")
    writeXStringSet( seq$dna , fastaName )
    
    gffName = paste0( "isoform_ts",ts,"_td",td,".gff")
    names(exonList) = paste( names(exonList) , translated , sep="_" )
    writeGff( isoform=exonList , file = gffName )
    #writeXStringSet( seq$dna[seq$translated] , fastaName )
    
    ###############################################################################################################
    kallisto = kallistoQuant( fastaName , SR1 , SR2 , EoutputDir )
    
    Sexp = log10(kallisto$tpm+1)
    Lexp = log10(correction[['chr2']]$normalizedIsoformCount+1)
    LSQuantCorr = cor.test(Lexp,Sexp)$estimate
    LSQuantPval = cor.test(Lexp,Sexp)$p.val
    
    ###############################################################################################################
    pdf( paste0( "Quant_LR_ts",ts,"_td",td,".pdf") )
    cols = sapply( seq$translated , function(x) ifelse(x,2,1) )
    plot( Lexp , Sexp , pch=16 , col=cols , main=paste0("Quantification by Long and Short Reads: r=",signif(LSQuantCorr,3)) ,  xlab="Log10 Long Read" , ylab="Log10 Short Read"  )
    abline(lm( Sexp~Lexp ))
    dev.off()
    
    ###############################################################################################################
    pdf( "gridSeach.pdf" )
    heatmap( score , Rowv = NA, Colv = NA, scale='none' )
    dev.off()
    
    ###############################################################################################################
    isoformNum = sum(sapply(correction,function(x)x$uniNum))
    isoformFrac = mean(sapply(correction,function(x)x$frac))
    info = data.frame( shortRead=ts , distance=td , isoformNum=isoformNum , isoformFrac=isoformFrac , translated=sum(seq$translated) , juncCorr , LSQuantCorr , LSQuantPval )
    write.table(info,"summary.txt",quote=F,sep="\t",col.names=T,row.names=F)
    
    
}


genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
genomeFasta = "/hpc/users/xzhus01/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"

# parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"

adjustNCjunc = TRUE
fixedMatchedLS = FALSE
useSJout = FALSE #TRUE 
fuzzyMatch = 100

chrom = "chr2"
s = 50147488
e = 51259537
cores = 30
thresSR=c(1:100) 
thresDis=c(1:30)
fixedMatchedLS=FALSE
LRflnc=NULL
LRnfl=NULL
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
system( paste("mkdir -p",folder) )


########################################################################################################################
##########################################   hipsc and adult  ####################################################################
########################################################################################################################

# 553 GABA
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027949/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/1_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/1_R2_001.fastq.gz"
outputDir=paste0(folder,"553_GABA")
SoutputDir = paste0(outputDir,"/SR")
LoutputDir = paste0(outputDir,"/LR")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores ,
                 SoutputDir=SoutputDir , LoutputDir=LoutputDir )


# 553 NGN2
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027950/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/2_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/2_R2_001.fastq.gz"
outputDir=paste0(folder,"553_NGN2")
SoutputDir = paste0(outputDir,"/SR")
LoutputDir = paste0(outputDir,"/LR")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores ,
                 SoutputDir=SoutputDir , LoutputDir=LoutputDir )


# fetal 23 weeks
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027966/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/6_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/6_R2_001.fastq.gz"
outputDir=paste0(folder,"fetal_23wks")
SoutputDir = paste0(outputDir,"/SR")
LoutputDir = paste0(outputDir,"/LR")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores ,
                 SoutputDir=SoutputDir , LoutputDir=LoutputDir )


# fetal 3 weeks
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027957/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/7_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/7_R2_001.fastq.gz"
outputDir=paste0(folder,"fetal_3wks")
SoutputDir = paste0(outputDir,"/SR")
LoutputDir = paste0(outputDir,"/LR")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores ,
                 SoutputDir=SoutputDir , LoutputDir=LoutputDir )


# 641 NGN2 
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027959/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/12_R1_001.fastq.gz"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/GENEWIZ/12_R2_001.fastq.gz"
outputDir = paste0(folder,"641_NGN2")
SoutputDir = paste0(outputDir,"/SR")
LoutputDir = paste0(outputDir,"/LR")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores ,
                 SoutputDir=SoutputDir , LoutputDir=LoutputDir )

# 641 NGN2 tofu
fasta="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027959/polished_high_qv_consensus_isoforms.fasta"
sam="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/641_NGN2/LR/Aligned.out.sam"
outputDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/641_NGN2/LR"
starTofu $fasta $sam $outputDir


# 581 NGN2
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027963/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R1.fastq"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R2.fastq"
outputDir = paste0(folder,"581_NGN2_fbnSR")
SoutputDir = paste0(outputDir,"/SR")
LoutputDir = paste0(outputDir,"/LR")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores ,
                 SoutputDir=SoutputDir , LoutputDir=LoutputDir )

# 581 NGN2 tofu
fasta="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027963/polished_high_qv_consensus_isoforms.fasta"
sam="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/581_NGN2_fbnSR/LR/Aligned.out.sam"
outputDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/581_NGN2_fbnSR/LR"
starTofu $fasta $sam $outputDir



# 581 GABA low number of phqv
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027964/polished_high_qv_consensus_isoforms.fasta"
SR1="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R1.fastq"
SR2="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R2.fastq"
outputDir = paste0(folder,"581_GABA_fbnSR")
SoutputDir = paste0(outputDir,"/SR")
LoutputDir = paste0(outputDir,"/LR")
STAR2bSMRT_NRXN( genomeDir=genomeDir, genomeFasta=genomeFasta, LRphqv=LRphqv, LRflnc=NULL, LRnfl=NULL,
                 SR1=SR1, SR2=SR2, useSJout=useSJout,  adjustNCjunc=adjustNCjunc, 
                 thresSR=thresSR, thresDis=thresDis, outputDir=outputDir, 
                 fixedMatchedLS=fixedMatchedLS, fuzzyMatch=fuzzyMatch, 
                 chrom=chrom , s=s , e=e , cores=cores ,
                 SoutputDir=SoutputDir , LoutputDir=LoutputDir )

# 581 GABA tofu low number of phqv
fasta="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027964/polished_high_qv_consensus_isoforms.fasta"
sam="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/581_GABA_fbnSR/LR/Aligned.out.sam"
outputDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/581_GABA_fbnSR/LR"
starTofu $fasta $sam $outputDir




# 641 GABA
LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/028024/polished_high_qv_consensus_isoforms.fasta"
outputDir = paste0(folder,"641_GABA")
SoutputDir = paste0(outputDir,"/SR")
LoutputDir = paste0(outputDir,"/LR")
starLong( genomeDir=genomeDir , LR=LRphqv , outputDir=LoutputDir , cores=20 , SJ=NULL )

# 641 GABA tofu
fasta="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/028024/polished_high_qv_consensus_isoforms.fasta"
sam="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/641_GABA/LR/Aligned.out.sam"
outputDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/641_GABA/LR"
starTofu $fasta $sam $outputDir


