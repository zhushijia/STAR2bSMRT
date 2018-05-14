library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
library(Biostrings)
library(foreach)
library(doMC)

cores = 30
thresSR=c(1:100) 
thresDis=c(1:30)
adjustNCjunc=FALSE
fixedMatchedLS=FALSE
chrom="chr17"
s= 90036900
e = 91089605
registerDoMC(cores)


outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/SRR1184043"

LoutputDir = paste0(outputDir,"/LR")
SoutputDir = paste0(outputDir,"/SR")
EoutputDir = paste0(outputDir,"/Exp")

system( paste0( "mkdir -p " , LoutputDir ) )
system( paste0( "mkdir -p " , SoutputDir ) )
system( paste0( "mkdir -p " , EoutputDir ) )

SRalignment = paste0(SoutputDir,"/alignments.bam")
LRalignment = paste0(LoutputDir,"/SRR1184043.sam")


SRjunc = getJunc( SRalignment , SoutputDir , chrom="chr17" , s= 90036900 , e = 91089605 )
LRinfo = getLRinfo3( LRalignment , NULL , LoutputDir , chrom="chr17" , s= 90036900 , e = 91089605 , jI=FALSE)
LRread = LRinfo$LRread
LRjunc = LRinfo$LRjunc
LRtag = LRinfo$LRtag

if( fixedMatchedLS )
{
  matchedLS = matchLSjunc( LRjunc , SRjunc )
} else {
  matchedLS = NULL
}

score = gridSearch( LRjunc , SRjunc , thresSR , thresDis , adjustNCjunc , matchedLS )

ij = which( score==max(score) , arr.ind=T )
ts = thresSR[ ij[1,1] ]
td = thresDis[ ij[1,2] ]
cat( ts , td , score[ij] , '\n ')

correction = generateCorrectedIsoform( LRjunc , SRjunc, LRtag , LRread  , ts , td , matchedLS )

setwd( EoutputDir )

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
gffName = paste0( "isoform_ts",ts,"_td",td,".gff")
exonList = juncToExon( juncList=correction[[chrom]]$isoform , s=50149082 , e=51255411 , exp=correction[[chrom]]$exp )
#writeGff( isoform=correction[[chrom]]$isoform , file = gffName , exp=correction[[chrom]]$exp , chrom='chr2' , s=50149082 , e=51255411 )
writeGff( isoform=exonList , file = gffName )

###############################################################################################################
#seq = generateSeq( genome=genome , isoform=correction[[chrom]]$isoform , exp=correction[[chrom]]$exp , chrom='chr2' , s=50149082 , e=51255411  )
seq = generateSeq( genome=genome , isoform=exonList )
fastaName = paste0( "isoform_ts",ts,"_td",td,".fa")
writeXStringSet( seq$dna , fastaName )
#writeXStringSet( seq$dna[seq$translated] , fastaName )

###############################################################################################################
kallisto = kallistoQuant( fastaName , SR1 , SR2 , EoutputDir )

Sexp = log10(kallisto$tpm+1)
Lexp = log10(correction[['chr2']]$exp+1)
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
isoformNum = sum(sapply(correction,function(x)x$num))
isoformFrac = mean(sapply(correction,function(x)x$frac))
info = data.frame( shortRead=ts , distance=td , isoformNum=isoformNum , isoformFrac=isoformFrac , translated=sum(seq$translated) , juncCorr , LSQuantCorr , LSQuantPval )
write.table(info,"summary.txt",quote=F,sep="\t",col.names=T,row.names=F)


