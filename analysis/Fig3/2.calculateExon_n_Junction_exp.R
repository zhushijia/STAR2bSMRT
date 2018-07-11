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
  rownames(fracs) = NULL
  fracs = data.frame(fracs)
  
  fracs$exon3b[ fracs$exon3a==1 & fracs$exon3b==1 ] = 0
  fracs$exon7a[ fracs$exon7a==1 & fracs$exon7b==1 ] = 0
  fracs$exon23a[ fracs$exon23a==1 & fracs$exon23b==1 ] = 0
  
  fracs
  
}

calculateExonExp = function(cc_fracs,ccExp)
{
  exonExpExpEachTranscript = apply( cc_fracs , 2 , function(x) x*rowSums(ccExp) )
  colSums(exonExpExpEachTranscript)
}

calculateJuncFrac = function(cc_fracs,ccExp)
{
  isoform = apply(cc_fracs , 1 , function(x) {
    iso1 = paste( which(x>0), collapse = "_" )
    paste0( "_" , iso1 , "_")
  })
  
  
  junc = matrix(ncol=ncol(cc_fracs),nrow=ncol(cc_fracs),data=0)
  colnames(junc) = colnames(cc_fracs)
  rownames(junc) = colnames(cc_fracs)
  
  for(i in 1:ncol(cc_fracs))
  {
    for(j in i:ncol(cc_fracs))
    {
      if(i==j)
      {
        junc[i,i] = sum( grepl( paste("",i,"",sep="_") , isoform ) )
      } else {
        junc[i,j] = sum( grepl( paste("",i,j,sep="_") , isoform ) )
      }
    }
    
  }
  
  frac = junc/diag(junc)
  frac[lower.tri(frac)] = -1
  diag(frac) = 0
  frac
}







annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)

library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
samples= c("581","641","2607","553")#,"642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal")
sampleNames = c("3' Del 1","3' Del 2","Control 1","Control 2")#,"Control 3","Adult 1","Adult 2","Adult 3","Fetal")

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
  
  gffs = gffs[1:4]
  
  translated = lapply( gffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
  }   )
  
  translatedGffs = lapply( 1:length(gffs) , function(i) gffs[[i]][ translated[[i]] ]  )
  names(translatedGffs) = names(gffs)
  
  translatedExps = lapply( translatedGffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
  }   )
  names(translatedExps) = names(gffs)
  
  
  caseGff = unionGff(translatedGffs[[1]],translatedGffs[[2]])
  contGff = unionGff(translatedGffs[[3]],translatedGffs[[4]])
  case_fracs = getFrac(caseGff,annotation)
  cont_fracs = getFrac(contGff,annotation)
  
  case_fracs2 = as.data.frame(round(case_fracs))
  cont_fracs2 = as.data.frame(round(cont_fracs))
  
  caseExp = sapply( 1:2 , function(i) {
    ee = translatedExps[[i]][ matchGff( caseGff , translatedGffs[[i]] ) ]
    ee[is.na(ee)] = 0
    ee
  } )
  contExp = sapply( 3:4 , function(i) {
    ee = translatedExps[[i]][ matchGff( contGff , translatedGffs[[i]] ) ]
    ee[is.na(ee)] = 0
    ee
  } )
  
  
  caseExonExp = calculateExonExp(case_fracs,caseExp)
  contExonExp = calculateExonExp(cont_fracs,contExp)
  
  caseJuncFrac = calculateJuncFrac(case_fracs,caseExp)
  contJuncFrac = calculateJuncFrac(cont_fracs,contExp)
  
  
  outputDir = paste0(folder,"/comparisonResult","/",parai)
  system( paste( "mkdir -p" , outputDir ) )
  setwd(outputDir)
  
  library(gplots)
  pdf("case.exon_junc_exp.pdf")
  barplot(caseExonExp,border=F,las=3)
  barplot( log10(caseExonExp) ,border=F,las=3)
  color = list( redgreen(256) , topo.colors(256) , cm.colors(256) , terrain.colors(256) , rainbow(256), heat.colors(256)  )
  heatmap.2(caseJuncFrac, col=color[[3]], scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
  dev.off()
  
  
  pdf("cont.exon_junc_exp.pdf")
  barplot(contExonExp,border=F,las=3)
  barplot( log10(contExonExp),border=F,las=3)
  color = list( redgreen(256) , topo.colors(256) , cm.colors(256) , terrain.colors(256) , rainbow(256), heat.colors(256)  )
  heatmap.2(contJuncFrac, col=color[[3]], scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
  dev.off()
  
  

}

















