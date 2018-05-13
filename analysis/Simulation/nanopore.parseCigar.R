library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
library(Biostrings)
library(foreach)
library(doMC)
library(data.table)

setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/minimap2/SRR6058583_NatMethod_E2")
sam <- file("SRR6058583_NatMethod_E2.sam", "r", blocking = FALSE)
#sam = fread("SRR6058583_NatMethod_E2.sam",quote="")
thres = 20

num=0
read = readLines(sam,n=1)
while( length(read)>0 )
{
  num = num+1
  cat(read,"\n")
  sep = strsplit(read,"\t")[[1]]
  if(length(sep)<10)
  {
    read = readLines(sam,n=1)
  } else {
    break
  }
}


juncs = list()

while( length(read)>0 )
{
  sep = strsplit(read,"\t")[[1]]
  start = as.integer(sep[4])
  cigar = parseCigar(sep[6])
  cigar = subset( cigar , !Op%in%c("I","S") )
  index = which( cigar[,1]>thres & as.character(cigar[,2]) %in% c("N","D") )
  
  juncL = c()
  juncR = c()
  
  for(i in 1:length(index))
  {
    juncL[i] = start + sum(cigar$num[1:(index[i]-1)])
    juncR[i] = start + sum(cigar$num[1:index[i]]) - 1
  }
  
  junc = paste(apply(data.frame(juncL,juncR),1,function(x)paste(x,collapse=",")) , collapse="," )
  juncs = c(juncs,junc)
  
  if(0)
  {
    junc1 = gsub( "^jI:B:i,","",sep[19] )
    junc2 = paste(apply(data.frame(juncL,juncR),1,function(x)paste(x,collapse=",")) , collapse="," )
    
    num = num+1
    cat(num,junc1==junc2,"\n")
    
    if( junc1!=junc2 )
      break
  }
  read = readLines(sam,n=1)
  read
}

close(sam)


