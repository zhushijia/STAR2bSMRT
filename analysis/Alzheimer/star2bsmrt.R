library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
library(Biostrings)
library(foreach)
library(doMC)

cores = 30
thresSR=c(1:10) 
thresDis=c(1:100)
adjustNCjunc=TRUE  
fixedMatchedLS=FALSE
thres = 20
registerDoMC(cores)


#LoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/intermediate_files/flnc_starlongNew"
LoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/intermediate_files/flnc_starlongNew/separateSam"
SoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/BioChain_A703252_Illumina/starShort"
EoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/Exp"

system( paste0( "mkdir -p " , LoutputDir ) )
system( paste0( "mkdir -p " , SoutputDir ) )
system( paste0( "mkdir -p " , EoutputDir ) )

SRalignment = paste0(SoutputDir,"/alignments.bam")

#LRalignment = paste0(LoutputDir,"/Aligned.out.sam")
LRalignment = paste0(LoutputDir,"/m150121_234021_42142_c100753522550000001823145207091506_s1_p0.sam")

#SRjunc = getJunc( SRalignment , SoutputDir  )
SJ.out.tab = read.table( paste0(SoutputDir,"/SJ.out.tab") , sep="\t")
SJ.out.tab = with( SJ.out.tab , data.frame(count=V7, chr=V1, start=V2, end=V3, motif=V5) )
SRjunc = split(SJ.out.tab , as.character(SJ.out.tab$chr) )


LRinfo = getLRinfo3( LRalignment , NULL , LoutputDir  , jI=TRUE)
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
sapply(correction,function(x)x$frac)
sapply(correction,function(x)x$num)

i = i+1
x0 = correction[[i]]$LSjuncCount
x1 = gff[grepl(paste0("SIRV",i),names(gff))]
x1 = sort(unique(do.call(c,lapply( x1 , function(x) if(nrow(x)>1) paste( x[1:(nrow(x)-1),3]+1 , x[2:nrow(x),2]-1 ,sep=",") ))))
x2 = sort(unique(do.call(c,lapply(correction[[i]]$isoform,function(x)paste(x[,2],x[,3],sep=",")))))
x1
x2
setdiff(x2,x1)




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
