library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

library(Biostrings)
library(foreach)
library(doMC)


genomeFasta = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
chrom = "chr2"
s = 50147488
e = 51259537
cores = 30
registerDoMC(cores)

thresSR=c(1:100) 
thresDis=c(1:30)
adjustNCjunc=FALSE
fixedMatchedLS=TRUE
folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/"
system( paste("mkdir -p",folder) )
genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
outputDir=paste0(folder,"553")

LoutputDir = paste0(outputDir,"/LR")
SoutputDir = paste0(outputDir,"/SR")
EoutputDir = paste0(outputDir,"/Exp")

system( paste0( "mkdir -p " , LoutputDir ) )
system( paste0( "mkdir -p " , SoutputDir ) )
system( paste0( "mkdir -p " , EoutputDir ) )

SRalignment = paste0(SoutputDir,"/alignments.bam")
LRalignment = paste0(LoutputDir,"/Aligned.out.sam")

LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24459_553/polished_high_qv_consensus_isoforms.fasta"
SRjunc = getJunc( SRalignment , SoutputDir , chrom , s= 50147488 , e = 51259537 )
LRinfo = getLRinfo( LRalignment , LR , LoutputDir , chrom , s= 50147488 , e = 51259537 )
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
LSjuncCount = subset(correction$chr2$LSjuncCount,lrCount>=50)
sj = LSjuncCount$srCount
lj = LSjuncCount$lrCount

fit = glm( sj ~ log(lj) , family=poisson() )
res = residuals(fit,type="response" )
predicts = predict(fit,type="response" )
data.frame(sj,lj,predicts,res)



fit = lm( log(sj) ~ log(lj) )
res = residuals(fit)
predicts = predict(fit) + min(res)
data.frame( sj , lj ,exp(predicts), as.integer(sj-exp(predicts)) )

data.frame(sj,lj,exp(predicts),res)










