library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
library(Biostrings)
library(foreach)
library(doMC)

cores = 30
thresSR=c(1:500) 
thresDis=c(1:30)
adjustNCjunc=TRUE  # TRUE and FALSE ts=408 td=15 num=281 frac=0.5069613
fixedMatchedLS=FALSE
chrom="chr17"
s= 90036900
e = 91089605
registerDoMC(cores)

annotation = read.table("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/Mouse_NRXN1_Annotation.txt",sep="\t",header=T)
annotation[6,2] = 90704367
annot_juncs = data.frame( js=annotation[2:nrow(annotation),2]+1 , je=annotation[1:(nrow(annotation)-1),1]-1 )

outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT"
LoutputDir = paste0(outputDir,"/LR")
SoutputDir = paste0(outputDir,"/SR")
EoutputDir = paste0(outputDir,"/Exp")

system( paste0( "mkdir -p " , LoutputDir ) )
system( paste0( "mkdir -p " , SoutputDir ) )
system( paste0( "mkdir -p " , EoutputDir ) )

SRalignment = paste0(SoutputDir,"/alignments.bam")
LRalignment = paste0(LoutputDir,"/Aligned.out.sam")


SRjunc = getJunc( SRalignment , SoutputDir , chrom="chr17" , s= 90036900 , e = 91089605 )
LRinfo = getLRinfo3( LRalignment , NULL , LoutputDir , chrom="chr17" , s= 90036900 , e = 91089605 , jI=TRUE)
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

juncs = do.call(c,lapply(correction[[1]]$isoform[correction[[1]]$exp>2],function(x)paste(x[,2],x[,3],sep=",")))
table(juncs)
sort(unique(juncs))


sites = c(sapply(strsplit(sort(unique(juncs)),","),function(x) as.integer(x)))


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

gffName = paste0( "isoform_ts",ts,"_td",td,".gff")
exonList = juncToExon( juncList=correction[[chrom]]$isoform , s=90036900 , e=91089605 , exp=correction[[chrom]]$exp )
writeGff( isoform=exonList , file = gffName )


#######################################################################################################
# check consistency with hiPSC

tag = lapply( exonList , function(x) apply(x,1,function(y) paste(y[2:3],collapse=",")) )
annot_tag = apply(annotation,1,function(y) paste(y[1:2],collapse=",") )

lapply(tag,function(x) as.character(annotation$Annotation)[ match(x,annot_tag) ]   )





