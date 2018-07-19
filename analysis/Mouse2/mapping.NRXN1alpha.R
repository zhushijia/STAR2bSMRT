startCondon = 90037077; endCondon = 91088726
startExon = 90036900; endExon = 91089605
startJunc =90037384; endJunc = 91087954
genomeDir="/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/gencode/GRCm38.p5/genome/StarIndexUsingPrimaryGtf"
LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/SRR1184043.fa"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/mouse-dIPFC-redo_S13_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/mouse-dIPFC-redo_S13_R2_001.fastq.gz"

library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")
library(Biostrings)
library(foreach)
library(doMC)

cores = 30
thresSR=c(1:100) 
thresDis=c(1:20)
adjustNCjunc=TRUE  
fixedMatchedLS=FALSE
fuzzyMatch=100
chrom="chr17"
s= 90036900
e = 91089605
registerDoMC(cores)

# adjustNCjunc=TRUE; fixedMatchedLS=FALSE; fuzzyMatch=100; getJuncBySJout
# ts=70 td=13 cor=0.9662687 uniNum=294 frac=0.7953995 for 3304 full-length reads



outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT"
LoutputDir = paste0(outputDir,"/LR_SRR1184043_STARlongNew")
SoutputDir = paste0(outputDir,"/SR")
EoutputDir = paste0(outputDir,"/Exp_SRR1184043_STARlongNew")

system( paste0( "mkdir -p " , LoutputDir ) )
system( paste0( "mkdir -p " , SoutputDir ) )
system( paste0( "mkdir -p " , EoutputDir ) )

SRalignment = paste0(SoutputDir,"/alignments.bam")
LRalignment = paste0(LoutputDir,"/Aligned.out.sam")


#SRjunc = getJuncBySam( SRalignment , SoutputDir , chrom="chr17" , s= 90036900 , e = 91089605 )
SRjunc = getJuncBySJout( SJout="SJ.out.tab" , SoutputDir , chrom="chr17" , s= 90036900 , e = 91089605 )

LRread = getReadByJI( LRalignment, LoutputDir )
LRread = subset( LRread , chr=='chr17' & start<90037384 & end>91087954 )
LRinfo = getLRinfo( LRread,  chrom=chrom, s=s, e=e )
LRread = LRinfo$LRread
LRjunc = LRinfo$LRjunc
LRtag = LRinfo$LRtag




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


#######################################################################################################
setwd( EoutputDir )
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


#######################################################################################################


tag = paste0( "exp", correction[[chrom]]$normalizedIsoformCount )
exonList = juncToExon( juncList=correction[[chrom]]$isoform , s=90037077 , e=91088726 , tag=tag )
genomeFasta = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/gencode/GRCm38.p5/genome/GRCm38.primary_assembly.genome.fa"
genome = readDNAStringSet(genomeFasta)
names(genome) = sapply( strsplit( names(genome) , ' ' ) , function(x) x[1] )
seq = generateSeq( genome=genome , isoform=exonList )
fastaName = paste0( "isoform_ts",ts,"_td",td,".fa")
translated = sapply( seq$translated , function(x) ifelse( x , "translated" , "untranslated" )  )
names(seq$dna) = paste(names(seq$dna),translated,sep="_")
writeXStringSet( seq$dna , fastaName )

gffName = paste0( "isoform_ts",ts,"_td",td,".gff")
names(exonList) = paste( names(exonList) , translated , sep="_" )
writeGff( isoform=exonList , file = gffName )



