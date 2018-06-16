library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone3/setup")
library(Biostrings)
library(foreach)
library(doMC)


setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a")
genome = readDNAStringSet("SIRV_isoforms_multi-fasta_170612a.fasta")
gff = readGff("SIRV_isoforms_multi-fasta-annotation_C_170612a.gff",chrom=NULL,s=0,e=Inf)

cores = 30
thresSR=c(1:10) 
thresDis=c(1:100)
adjustNCjunc=TRUE  # TRUE and FALSE ts=408 td=15 num=281 frac=0.5069613
fixedMatchedLS=FALSE
thres = 20
registerDoMC(cores)


#LoutputDir = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/minimap2/SRR6058583_NatMethod_E2"
LoutputDir = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/minimap2/SRR5286959_1"
SoutputDir = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/starLong/SRR5286956_Lexogen_100fg_2_Smartseq2_ONT/SR"
EoutputDir = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/starLong/SRR5286956_Lexogen_100fg_2_Smartseq2_ONT/Exp"

system( paste0( "mkdir -p " , LoutputDir ) )
system( paste0( "mkdir -p " , SoutputDir ) )
system( paste0( "mkdir -p " , EoutputDir ) )

SRalignment = paste0(SoutputDir,"/alignments.bam")
#LRalignment = paste0(LoutputDir,"/SRR6058583_NatMethod_E2.sam")
LRalignment = paste0(LoutputDir,"/SRR5286959.sam")


SRjunc = getJuncBySam( SRalignment , SoutputDir  )
#SRjunc = getJuncBySJout( "SJ.out.tab" , SoutputDir  )

LRinfo = getLRinfo( LRalignment , NULL , LoutputDir  , jI=FALSE)
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
setdiff(x2,x1) # false positive
setdiff(x1,x2) # false negative: what we miss




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
