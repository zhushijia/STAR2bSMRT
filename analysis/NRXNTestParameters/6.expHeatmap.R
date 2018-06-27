library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

#annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1.txt',sep='\t')
annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1ExonAnnotations.txt',sep='\t',header=T)
annotation[2,2] = annotation[3,3]+1
annotation[25,3] = annotation[24,2]-1

######################################
annotation[27,] = annotation[8,]
annotation[8,3] = 50848363
annotation[27,2] = 50848364
annotation[,1] = as.character(annotation[,1])
annotation[27,1] = "exon7c"
annotation[27,4] = 23

library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone4/setup")

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
samples= c("2607","553","581","641","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal")
sampleNames = c("Control 1","Control 2","3' Del 1","3' Del 2","Control 3","Adult 1","Adult 2","Adult 3","Fetal")


parameters = dir( paste0(folder,"/fetal") )

for(parai in parameters)  
{
  
  cat(parai,"\n")
  
  gffs = lapply( samples , function(sample)
  {
  	files = dir( paste0(folder,"/",sample,"/",parai) )
  	gff = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
  	readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
  } )
  names(gffs) = sampleNames
  
  translated = lapply( gffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
  }   )
  
  
  translatedGffs = lapply( 1:length(gffs) , function(i) gffs[[i]][ translated[[i]] ]  )
  names(translatedGffs) = sampleNames
  
  hiPSC = unionGff(translatedGffs[[1]],translatedGffs[[2]])
  patient = unionGff(translatedGffs[[3]],translatedGffs[[4]])
  mergeTranslatedGffs = list(  case=patient ,control=hiPSC )
  mergeTranslatedExps = lapply( mergeTranslatedGffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
  }   )

  
  getFrac = function(gff,annotation)
  {
    
    frac = function(gs,annot_sites)
    {
      sapply(annot_sites,function(x) {
        max( sapply(gs,function(y) mean(x%in%y)) )
      } )
    }
    
    gff_sites = lapply(gff,function(y) apply(y,1,function(z) z[2]:z[3] ) )
    annot_sites = apply(annotation,1,function(z) z[2]:z[3] )
    fracs = do.call(rbind,lapply( gff_sites , function(gs) frac(gs,annot_sites) ) )
    colnames(fracs) = as.character(annotation$Exon)
    fracs
  }
  
  caseExp = mergeTranslatedExps[[1]]
  contExp = mergeTranslatedExps[[2]]
  caseGff = mergeTranslatedGffs[[1]]
  contGff = mergeTranslatedGffs[[2]]
  case_fracs2 = getFrac(caseGff,annotation)
  cont_fracs2 = getFrac(contGff,annotation)
  range = matchGff(contGff,caseGff)
  range = range[!is.na(range)]
  case_fracs2[range, ] = -case_fracs2[range, ]
  range = matchGff(caseGff,contGff)
  range = range[!is.na(range)]
  cont_fracs2[range, ] = -cont_fracs2[range, ]
  
  plotHeatmap = function(cc_fracs,ccExp)
  {
    require('gplots')
    ccExp = log(ccExp+1)
    x = cc_fracs
    x = cbind( (cc_fracs) , 0 , exp=1-ccExp/max(ccExp) )
    x = x[order(x[,29]),]
    #heatmap( x ,Rowv=NA,Colv=NA , scale="none" , col=redgreen(256) )
    heatmap.2(x, col=cm.colors(256), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
  }
  
  outputDir = paste0(folder,"/comparisonResult","/",parai)
  system( paste( "mkdir -p" , outputDir ) )
  setwd(outputDir)
  
  pdf("case.heatmap_translated_ConsiderOverlap.pdf")
  plotHeatmap( case_fracs2 , caseExp )
  dev.off()
  
  pdf("cont.heatmap_translated_ConsiderOverlap.pdf")
  plotHeatmap( cont_fracs2 , contExp )
  dev.off()
  

  
}
