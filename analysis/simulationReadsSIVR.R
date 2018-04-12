library(polyester)
library(Biostrings)
library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a")
outputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/simulation/E2_simulation"
fasta_file = "SIRV_isoforms_multi-fasta-annotation_C_170612a_shijia.fa"
readlen = 50

fasta = readDNAStringSet(fasta_file)
exp = 
c( c( 1, 1/4, 1, 4, 1, 4, 1/4, 1/32),
c( 1/4, 4, 1/4, 1, 1/32, 1/4),
c( 4, 1/32, 1, 1/32, 1/4, 1/4, 1, 4, 1/4, 1/32, 4),
c( 1/32, 1, 1/32, 4, 1/32, 1/4, 4),
c( 1/4, 1/4, 1/4, 1, 1, 1/32, 4, 1/4, 1/32, 1/32, 1/4), # no 502: 1/32
c( 1/32, 4, 1/32, 1, 1/4, 1/4, 1, 1/4, 1, 1/4, 1/32, 4, 1, 1/4, 1/4, 4, 1, 1/32), 
c( 1/4, 1/32, 1/32, 1/32, 1, 1, 1 ) )

fold_changes = matrix( rep(1,2*length(fasta)) , nrow=length(fasta) )
setwd(outputDir)
writeXStringSet( fasta, 'E2.fa')
readspertx = round( 100*exp * width(fasta) / readlen )

# simulation call:
simulate_experiment('E2.fa', readlen=readlen , reads_per_transcript=readspertx, 
                    num_reps=c(1,1), fold_changes=fold_changes, outdir='simulated_reads') 

gff=readGff("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a/SIRV_isoforms_multi-fasta-annotation_C_170612a.gff",chrom=NULL,s=0,e=Inf)
junc = lapply( 1:length(gff), function(i) {
  x = gff[[i]]
  e = exp[i]
  tag = paste( x[-nrow(x),1] , x[-nrow(x),3]+1 , x[-1,2]-1 )
  if(length(tag)>0)
  {
    data.frame(tag,e)
  } else {
    NULL
  }
} )
LRjunc1 = do.call(rbind,junc)
x = tapply(LRjunc1[,2],LRjunc1[,1],sum)
LRjunc1 = data.frame( tag=names(x) , e=x )
LRjunc2 = do.call(rbind,LRjunc)
LRjunc2 = data.frame( tag=paste(LRjunc2[,2],LRjunc2[,3],LRjunc2[,4]) , e=LRjunc2[,1] )
LRjunc2 = LRjunc2[ match( as.character(LRjunc1[,1]) , as.character(LRjunc2[,1]) ) ,  ]
LRjunc2 = LRjunc2[!is.na(LRjunc2$e),]
LRjunc1 = LRjunc1[ match( as.character(LRjunc2[,1]) , as.character(LRjunc1[,1]) ) ,  ]
cor( LRjunc1$e , LRjunc2$e )





library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
library(foreach)
library(doMC)
registerDoMC(20)
genomeFasta = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
cores=30
thresSR=c(1:100)
thresDis=c(1:30)
adjustNCjunc=TRUE
fixedMatchedLS=FALSE
genomeDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a/STAR_index"
SR1="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/simulation/E2_simulation/simulated_reads/sample_01_1.fasta"
SR2="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/simulation/E2_simulation/simulated_reads/sample_01_2.fasta"
outputDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/starLong/SRR5286959"

LoutputDir = paste0(outputDir,"/LR")
SoutputDir = paste0(outputDir,"/SR")
EoutputDir = paste0(outputDir,"/Exp")

system( paste0( "mkdir -p " , LoutputDir ) )
system( paste0( "mkdir -p " , SoutputDir ) )
system( paste0( "mkdir -p " , EoutputDir ) )

starShort( genomeDir , SR1 , SR2 , SoutputDir )

SRalignment = paste0(SoutputDir,"/alignments.bam")
LRalignment = paste0(LoutputDir,"/Aligned.out.sam")

SRjunc = getJunc( SRalignment , SoutputDir , chrom=NULL , s= 0 , e = Inf )
LRinfo = getLRinfo2( LRalignment , NULL , LoutputDir , chrom=NULL , s= 0 , e = Inf )
LRread = LRinfo$LRread
LRjunc = LRinfo$LRjunc
LRtag = LRinfo$LRtag

if( fixedMatchedLS )
{
  matchedLS = matchLSjunc( LRjunc , SRjunc )
} else {
  matchedLS = NULL
}

thresSR=c(0:10)
thresDis=c(0:100)
score = gridSearch( LRjunc , SRjunc , thresSR , thresDis , adjustNCjunc , matchedLS )

ij = which( score==max(score) , arr.ind=T )
ts = thresSR[ ij[1,1] ]
td = thresDis[ ij[1,2] ]
cat( ts , td , score[ij] , '\n ')


for (i in 1:7)
{
	gg = gff[chr==paste0("gene_id SIRV",i)]
	jj = sapply(gg,function(g) {
		paste( g[-nrow(g),3]+1 , g[-1,2]-1 )
		} )
	table(do.call(c,jj))
	subset(x,chr==paste0("SIRV",i))
	jj

}


