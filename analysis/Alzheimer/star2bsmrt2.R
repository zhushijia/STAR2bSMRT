library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone3/setup")
library(Biostrings)
library(foreach)
library(doMC)

cores = 30
thresSR=c(1:100) 
thresDis=c(1:30)
adjustNCjunc=TRUE  
fixedMatchedLS=FALSE
thres = 20
registerDoMC(cores)


FLNCoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/intermediate_files/flnc_starlongNew/"
NFLoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/intermediate_files/nfl_starlongNew/"

SoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/BioChain_A703252_Illumina/mapping_hg19"

system( paste0( "mkdir -p " , FLNCoutputDir ) )
system( paste0( "mkdir -p " , NFLoutputDir ) )
#system( paste0( "mkdir -p " , EoutputDir ) )

SRalignment = paste0(SoutputDir,"/alignments.bam")
FLNCalignment = paste0(FLNCoutputDir,"/Aligned.out.sam")
NFLalignment = paste0(NFLoutputDir,"/Aligned.out.sam")

SRjunc = getJuncBySJout( SJout="SJ.out.tab" , SoutputDir  )

phqv=NULL 
chrom=NULL
s=0
e=Inf
jI=TRUE
matchedLS = NULL
alignments=FLNCalignment 
outputDir=FLNCoutputDir
lr1 = getReadByJI( alignments , outputDir )

alignments=NFLalignment 
outputDir=NFLoutputDir
lr2 = getReadByJI( alignments , outputDir )

LRread = rbind(lr1,lr2) 

LRinfo = getLRinfo( alignment , phqv=NULL , outputDir , jI=TRUE)
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



zz = coef(nnnpls(as.matrix(x),y,con=rep(1,ncol(x))))
fit = as.matrix(x)%*%as.matrix(zz)
cor(fit,y)


#splitSmrtcell( alignments=LRalignment , outputDir=LoutputDir )


LRinfo = list()
for( sc in dir( paste0(LoutputDir,"/separateSam/") )  )
{
  outputDir = paste0(LoutputDir,"/separateSam/",sc)
  alignment = paste0( outputDir,"/Aligned.out.sam")
  LRinfo[[sc]] = getLRinfo( alignment , phqv=NULL , outputDir , jI=TRUE)
}


LRjuncAll = lapply( LRinfo , function(x) do.call(rbind,x$LRjunc) )
LRjuncAll2 = Reduce( function(x, y) merge(x, y, by=c("chr","start","end"), all=TRUE) , LRjuncAll )
LRjuncAll2[is.na(LRjuncAll2)] = 0
LRjuncAll3 = split(LRjuncAll2,LRjuncAll2$chr)

#SRjunc = getJuncBySam( SRalignment , SoutputDir  )
SRjunc = getJuncBySJout( SJout="SJ.out.tab" , SoutputDir  )



ts=5
td=30
CHR = intersect( names(SRjunc) , names(LRjuncAll3) )
CHR = CHR[ grepl('chr',CHR)  ]
LSjuncCount = foreach( k=1:length(CHR) ) %dopar%
  
LSjuncCount = list()
for( k in 1:length(CHR) ) 
{
  chr = CHR[k]
  cat(chr,"\n")
  lrc = LRjuncAll3[[chr]]
  src = SRjunc[[chr]]
  src = subset( src , count>=ts )
  SRmatch = matchLSjuncOneChr( lrc , src )
  range = which( SRmatch[,'LSdis']<=td ) 
  LRcorres = data.frame(lrc,SRmatch)[range, ]
  lrCount = apply( LRcorres[,-c(1:3,ncol(LRcorres),ncol(LRcorres)-1) ] , 
                   2 , function(x) tapply( x , LRcorres$SRindex , sum ) )
  dat = data.frame( lrCount , src[as.integer(rownames(lrCount)),] )
  LSjuncCount[[chr]] = dat
}

dat = do.call(rbind,LSjuncCount)
x = dat[,1:40]
y = dat[,41]

z = lm( y~as.matrix(x)+1 ) 
cor( predict(z) , y  )
cor( predict(z) , y , method='spearman' )


z = lm( log(y+1)~log(as.matrix(x)+1) ) 
cor( predict(z) , log(y+1)  )
cor( predict(z) , log(y+1) , method='spearman' )


plot( predict(z) , log(y+1)  )
smoothScatter( log(predict(z)) , log(y+1)  )
dev.off()

a=as.list(rep(0,40))
names(a) = colnames(X)[1:40]
formula = as.formula( paste0( "X41~" , paste( colnames(X)[1:40] , collapse="+" ) ))
nls( formula , data=X , start=a , lower = 0, algorithm = "plinear")






z = lm( log(y+1)~log(as.matrix(x)+1) ) 
cor( predict(z) , log(y+1)  )
plot( predict(z) , log(y+1)  )
dev.off()

z <- nnls( as.matrix(x) , y )
coef(z)
fit <- as.vector(coef(z)%*%t(x))

nls( y~a*x1 + b*x2 +c*x3...,data=mydata, 
    start=c(a=1,b=1,c=1...), lower=c(a=0,b=0,c=NA,...), algorithm="port")




