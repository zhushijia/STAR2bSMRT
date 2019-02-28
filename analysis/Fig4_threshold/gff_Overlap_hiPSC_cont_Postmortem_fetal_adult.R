library(STAR2bSMRT,lib.loc="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone4/setup")

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters"
parai1 = "adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_TRUE_fuzzyMatch_100"
parai2 = "adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"

parai = parai2
samples = c( paste0( c("2607","553","642",
                       "fetal","fetal_23wks","fetal_3wks",
                       "adult_dlPFC1_10","NRXN_adult_dlPFC1_12","adult_dlPFC1_13"),"/",parai) )
sampleName = c("cont1","cont2","cont3",'fetal23w1','fetal23w2','fetal3w','adult1','adult2','adult3')

gff = lapply( samples , function(sample)
{
    files = dir( paste0(folder,"/",sample,"/") )
    tmp = paste0(folder,"/",sample,"/",files[grepl("gff$",files)])
    readGff( tmp , chrom='chr2' , s=50149082 , e=51255411 )
} )

exp = lapply( gff , function(x) {
    sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
}   )

gffs = lapply( 1:length(gff), function(i) gff[[i]][order(exp[[i]],decreasing=T)] )
exps = lapply( 1:length(gff), function(i) exp[[i]][order(exp[[i]],decreasing=T)] )
names(gffs) = sampleName
names(exps) = sampleName

translated = lapply( gffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
}   )


#############################################################################################
######### pca
#############################################################################################

fthres = function(thres)
{
  #thres = 8
  tGffs = lapply( 1:length(gffs) , function(i) gffs[[i]][ translated[[i]] & exps[[i]]>=thres ]  )
  tExps = lapply( 1:length(gffs) , function(i) exps[[i]][ translated[[i]] & exps[[i]]>=thres ]  )
  names(tGffs) = sampleName
  names(tExps) = sampleName
  
  #tGffs = tGffs[c(1:5,8:11)]
  #tExps = tExps[c(1:5,8:11)]
  sapply(tGffs,length)
  
  allGffs = tGffs[[1]]
  for(i in 2:length(tGffs))
    allGffs = unionGff(allGffs,tGffs[[i]])
  
  labels = names(tGffs)
  x1 = sapply( 1:length(tGffs) , function(i) {
    flag = rep(0,length(allGffs))
    index = matchGff( allGffs, tGffs[[i]] )
    flag[ !is.na(index) ] = 1
    flag
  })
  colnames(x1) = labels
  
  
  x2 = sapply( 1:length(tGffs) , function(i) {
    flag = rep(0,length(allGffs))
    index = matchGff( allGffs, tGffs[[i]] )
    flag[ !is.na(index) ] = tExps[[i]][ index[!is.na(index)] ]
    flag
  })
  colnames(x2) = labels
  
  
  pc = prcomp( t(x1) )$x
  plot( pc[,1] , pc[,2] , col="white",main=paste0("Isoform Yes or No: thres=",thres) )
  text( pc[,1] , pc[,2] , colnames(x1) )
  dd <- dist( t(x1) , method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  plot(hc,main=paste0("Isoform Yes or No: thres=",thres) )
  
  pc = prcomp( t(x2) )$x
  plot( pc[,1] , pc[,2] , col="white" ,main=paste0("Isoform Long Read Count: thres=",thres) )
  text( pc[,1] , pc[,2] , colnames(x2) )
  dd <- dist( t(x2) , method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  plot(hc,main=paste0("Isoform Long Read Count: thres=",thres) )
  

}


ftopn = function(topn)
{
  tGff = lapply( 1:length(gffs) , function(i) gffs[[i]][ translated[[i]] ]  )
  tExp = lapply( 1:length(gffs) , function(i) exps[[i]][ translated[[i]] ]  )
  
  tGffs = lapply( 1:length(gffs) , function(i) tGff[[i]][ 1:min(topn,length(tGff[[i]])) ]  )
  tExps = lapply( 1:length(gffs) , function(i) tExp[[i]][ 1:min(topn,length(tExp[[i]])) ]  )
  names(tGffs) = sampleName
  names(tExps) = sampleName
  
  #tGffs = tGffs[c(1:5,8:11)]
  #tExps = tExps[c(1:5,8:11)]
  sapply(tGffs,length)
  
  allGffs = tGffs[[1]]
  for(i in 2:length(tGffs))
    allGffs = unionGff(allGffs,tGffs[[i]])
  
  labels = names(tGffs)
  x1 = sapply( 1:length(tGffs) , function(i) {
    flag = rep(0,length(allGffs))
    index = matchGff( allGffs, tGffs[[i]] )
    flag[ !is.na(index) ] = 1
    flag
  })
  colnames(x1) = labels
  
  
  x2 = sapply( 1:length(tGffs) , function(i) {
    flag = rep(0,length(allGffs))
    index = matchGff( allGffs, tGffs[[i]] )
    flag[ !is.na(index) ] = tExps[[i]][ index[!is.na(index)] ]
    flag
  })
  colnames(x2) = labels
  
  
  pc = prcomp( t(x1) )$x
  plot( pc[,1] , pc[,2] , col="white",main=paste0("Isoform Yes or No: topn=",topn) )
  text( pc[,1] , pc[,2] , colnames(x1) )
  dd <- dist( t(x1) , method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  plot(hc,main=paste0("Isoform Yes or No: topn=",topn) )
  
  pc = prcomp( t(x2) )$x
  plot( pc[,1] , pc[,2] , col="white" ,main=paste0("Isoform Long Read Count: topn=",topn) )
  text( pc[,1] , pc[,2] , colnames(x2) )
  dd <- dist( t(x2) , method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  plot(hc,main=paste0("Isoform Long Read Count: topn=",topn) )
  
  
}


sapply(tGffs,length)


setwd( paste0(folder,"/comparisonResult/",parai ) )
pdf("PCA_new_cont_fetal_adult_hiPSC_topN.pdf")
for(topn in c(1:10)*10 )
{
  cat(topn,'\n')
  ftopn(topn)
}
dev.off()

pdf("PCA_new_cont_fetal_adult_hiPSC_thres.pdf")
for(thres in c(1:30) )
{
  cat(thres,'\n')
  fthres(thres)
}
dev.off()


