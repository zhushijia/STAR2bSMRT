library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
library(Biostrings)
library(foreach)
library(doMC)

cores = 30
thresSR=c(1:100) 
thresDis=c(1:30)
adjustNCjunc=TRUE
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
correction[[1]]$frac
correction[[1]]$num

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




t = t+1
x = tag[[t]]
y = x[ !x%in%correctTag ]
tt = paste(LRjunc[[1]][,3],LRjunc[[1]][,4],sep=",")
ind = SRmatch[tt%in%y , ]
ind
src[ind[,1],]
src[ind[1],]
y






