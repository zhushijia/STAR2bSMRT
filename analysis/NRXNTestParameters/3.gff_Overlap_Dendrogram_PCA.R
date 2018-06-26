library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
samples= dir()
samples= c("2607","553","581","641","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal")

parai = "adjustNCjunc_FALSE_fixedMatchedLS_FALSE_useSJout_TRUE_fuzzyMatch_0"
parai = "adjustNCjunc_FALSE_fixedMatchedLS_FALSE_useSJout_TRUE_fuzzyMatch_100"
parai = "adjustNCjunc_FALSE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_0"

gffs = lapply( samples , function(sample)
{
	files = dir( paste0(folder,"/",sample,"/",parai) )
	gff = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
	readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(gffs) = samples

exp = lapply( gffs , function(x) {
  sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
}   )

translated = lapply( gffs , function(x) {
  sapply( strsplit(names(x),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
}   )


translatedGffs = lapply( 1:length(gffs) , function(i) gffs[[i]][ translated[[i]] ]  )
names(translatedGffs) = samples

hiPSC = unionGff(translatedGffs[[1]],translatedGffs[[2]])
humanBrain = unionGff(unionGff(translatedGffs[[6]],translatedGffs[[7]]),translatedGffs[[8]])
fetal = translatedGffs[[9]]
patient = unionGff(translatedGffs[[3]],translatedGffs[[4]])
mergeTranslatedGffs = list( hiPSC=hiPSC , adult=humanBrain , fetal=fetal , patient=patient )


outputDir = paste0("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/comparisonResult","/",parai)
system( paste( "mkdir -p" , outputDir ) )
setwd(outputDir)

vennDiagramGff(translatedGffs[c(1:5)] , filename=paste0("vennDiagram_",5,".png") )
vennDiagramGff(translatedGffs[c(1:4)] , filename=paste0("vennDiagram_",4,".png") )
vennDiagramGff( mergeTranslatedGffs[1:3] , filename ="vennDiagram_hiPSC_humanbrain_translated.png" )


###########################################################################################
########################  dendrogram and pca ###################################
###########################################################################################

dendrogram_pca = function( gffs, fileName )
{

  allGffs = gffs[[1]]
  for(i in 2:length(gffs))
	allGffs = unionGff(allGffs,gffs[[i]])

	exp = lapply( gffs , function(x) {
	  sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
	}   )


  labels = names(gffs)
  x1 = sapply( 1:length(gffs) , function(i) {
    flag = rep(0,length(allGffs))
    ind = matchGff( allGffs, gffs[[i]] )
    flag[ !is.na(ind) ] = 1
    flag
  })
  colnames(x1) = labels
  
  x2 = sapply( 1:length(gffs) , function(i) {
    flag = rep(0,length(allGffs))
    ind = matchGff( allGffs, gffs[[i]] )
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

dendrogram_pca( translatedGffs[c(1:9)] , "dendrogram_translatedAllGffs.pdf" )
dendrogram_pca( translatedGffs[c(1:2,5:9)] , "dendrogram_translatedPartGffs.pdf" )

dendrogram_pca( gffs[c(1:9)] , "dendrogram_allGffs.pdf" )
dendrogram_pca( gffs[c(1:2,5:9)] , "dendrogram_partGffs.pdf" )

