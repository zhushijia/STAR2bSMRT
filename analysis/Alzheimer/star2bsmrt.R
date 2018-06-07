library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone3/setup")
library(Biostrings)
library(foreach)
library(doMC)

cores = 20
thresSR=c(1:100) 
thresDis=c(1:30)
adjustNCjunc=TRUE  
fixedMatchedLS=FALSE
thres = 20
registerDoMC(cores)


#LoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/intermediate_files/flnc_starlongNew"
#LoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/intermediate_files/flnc_starlongNew/separateSam/m141129_125831_42161_c100698142550000001823143403261592_s1_p0"
LoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/intermediate_files/flnc_starlongNew/"
SoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/BioChain_A703252_Illumina/mapping_hg19"
#EoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/intermediate_files/flnc_starlongNew/separateSam/m141129_125831_42161_c100698142550000001823143403261592_s1_p0/Exp"

system( paste0( "mkdir -p " , LoutputDir ) )
system( paste0( "mkdir -p " , SoutputDir ) )
#system( paste0( "mkdir -p " , EoutputDir ) )

SRalignment = paste0(SoutputDir,"/alignments.bam")
LRalignment = paste0(LoutputDir,"/Aligned.out.sam")

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



ts=2
td=15
chr="chr2"
CHR = intersect( names(SRjunc) , names(LRjuncAll3) )
CHR = CHR[ grepl('chr',CHR)  ]
LSjuncCount = foreach( k=1:length(CHR) ) %dopar%
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
  data.frame( lrCount , src[as.integer(rownames(lrCount)),] )
}


lrc = LRjuncAll3[[chr]]
src = SRjunc[[chr]]
src = subset( src , count>=ts )
SRmatch = matchLSjuncOneChr( lrc , src )
range = which( SRmatch[,'LSdis']<=td ) 
LRcorres = data.frame(lrc,SRmatch)[range, ]
lrCount = apply( LRcorres[,-c(1:3,ncol(LRcorres),ncol(LRcorres)-1) ] , 
                 2 , function(x) tapply( x , LRcorres$SRindex , sum ) )
dat = data.frame( lrCount , src[as.integer(rownames(lrCount)),] )
x = dat[,1:40]
y = dat[,41]


z = lm( y~as.matrix(x)+1 ) 
cor( predict(z) , y  )

z = lm( log(y+1)~log(as.matrix(x)+1) ) 
cor( predict(z) , log(y+1)  )


plot( predict(z) , log(y+1)  )
smoothScatter( predict(z) , log(y+1)  )
dev.off()









z = lm( log(y+1)~log(as.matrix(x)+1) ) 
cor( predict(z) , log(y+1)  )
plot( predict(z) , log(y+1)  )
dev.off()

z <- nnls( as.matrix(x) , y )
coef(z)
fit <- as.vector(coef(z)%*%t(x))

nls( y~a*x1 + b*x2 +c*x3...,data=mydata, 
    start=c(a=1,b=1,c=1...), lower=c(a=0,b=0,c=NA,...), algorithm="port")


z <- glm( log(y+1) ~ log(as.matrix(x)+1), family=poisson()) 
summary(z)
cor( predict(z) , y  )
cor( log(predict(z)+1) , log(y+1)  )




