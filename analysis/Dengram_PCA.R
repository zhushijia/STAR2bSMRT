library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/"
setwd(folder)
samples= c("2607","553","581","641","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal")

gffs = lapply( samples , function(sample)
{
  files = dir( paste0(sample,"/Exp") )
  gff = paste0(sample,"/Exp/",files[grep("gff",files)])
  readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(gffs) = samples

exp = lapply( gffs , function(x) {
  sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
}   )


allGff = gffs[[1]]
for(i in 2:length(gffs))
  allGff = unionGff(allGff,gffs[[i]])

library(Biostrings)
ref = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
genome = readDNAStringSet(ref)
allSeq = generateSeq( genome , isoform=allGff )


setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/result/pipeline_nonAdjustNCjunc")

dendrogram_pca = function( range , y , fileName )
{
  labels = names(gffs)[range]
  x1 = sapply( range , function(i) {
    flag = rep(0,length(y))
    ind = matchGff( y, gffs[[i]] )
    flag[ !is.na(ind) ] = 1
    flag
  })
  colnames(x1) = labels
  
  x2 = sapply( range , function(i) {
    flag = rep(0,length(y))
    ind = matchGff( y, gffs[[i]] )
    flag[ !is.na(ind) ] = exp[[i]][ ind[!is.na(ind)] ]
    flag
  })
  colnames(x2) = labels
  
  
  pdf(fileName)
  dd <- dist( t(x1) , method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  plot(hc,main="Isoform Yes or No")
  
  pc = prcomp( t(x1) )$x
  plot( pc[,1] , pc[,2] , col="white",main="Isoform Yes or No")
  text( pc[,1] , pc[,2] , colnames(x1) )
  
  dd <- dist( t(x2) , method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  plot(hc,main="Isoform Long Read Count")
  
  pc = prcomp( t(x2) )$x
  plot( pc[,1] , pc[,2] , col="white" ,main="Isoform Long Read Count")
  text( pc[,1] , pc[,2] , colnames(x2) )
  
  dev.off()
}

dendrogram_pca( c(1:9) , allGff , "dendrogram_allGffs.pdf" )
dendrogram_pca( c(1:9) , allGff , "dendrogram_allTranslatedGffs.pdf" )
dendrogram_pca( c(1:2,5:9) , allGff[allSeq$translated] , "dendrogram_partGffs.pdf" )
dendrogram_pca( c(1:2,5:9) , allGff[allSeq$translated] , "dendrogram_partTranslatedGffs.pdf" )




