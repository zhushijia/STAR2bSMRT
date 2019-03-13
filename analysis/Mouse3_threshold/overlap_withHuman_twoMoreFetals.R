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


fracToTag = function(fracs)
{
  apply(fracs,1,function(x) { 
    paste(colnames(fracs)[x>0.9] , collapse="_" )
  })
}


unionGffList = function(GffList)
{
  tmp = GffList[[1]]
  if( length(GffList) > 1 )
  {
    for (i in 2:length(GffList)) tmp = unionGff( tmp , GffList[[i]] )
  }
  tmp
}

getFracFromGffList = function(GffList,annotation)
{
  tmp = unionGffList(GffList)
  getFrac(tmp,annotation)
}

readHumanGff = function(samples,threshold)
{
  lapply( samples , function(sample)
  {
    files = dir( paste0(folder,"/",sample,"/",parai) )
    file = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
    gff = readGff( file , chrom='chr2' , s=50149082 , e=51255411 )
    translated = sapply( strsplit(names(gff),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
    gff = gff[translated]
    exp = sapply( strsplit(names(gff),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
    gff = gff[exp>=threshold]
    exp = exp[exp>=threshold]
    list(gff,exp)
  } )
}
  

#####################################################################################
##########                              ########
#####################################################################################

library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")
humanAnnotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
mouseAnnotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT/NRXN1_mm10_ExonAnnotations_shijia.txt',sep='\t',header=T)

#####################################################################################
##########       Human sample                      ########
#####################################################################################
threshold = 0

for(threshold in c(1:10))
{
  cat("threshold =",threshold,"\n")
  samples= c("581","641","2607","553","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks") #
  sampleNames = c("3' Del 1","3' Del 2","Control 1","Control 2","Adult 1","Adult 2","Adult 3","Fetal 1","Fetal 2","Fetal 3") #
  folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
  parai = "adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
  humanGffExp = readHumanGff(samples,threshold=threshold)
  names(humanGffExp) = sampleNames
  humanGffs = lapply(humanGffExp,function(x) x[[1]]  )
  humanExps = lapply(humanGffExp,function(x) x[[2]]  )
  
  #####################################################################################
  ##########      Mouse sample                     ########
  #####################################################################################
  
  file = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT/Exp_SRR1184043_STARlongNew/isoform_ts70_td13.gff"
  mouseGffs = readGff( file , chrom="chr17" , s= 90036900 , e = 91089605 )
  translated = sapply( strsplit(names(mouseGffs),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
  mouseGffs = mouseGffs[translated]
  mouseExps = sapply( strsplit(names(mouseGffs),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
  mouseGffs = mouseGffs[mouseExps>=threshold]
  mouseExps = mouseExps[mouseExps>=threshold]
  Gffs = c( humanGffs, Mouse=list(mouseGffs) )
  Exps = c( humanExps, Mouse=list(mouseExps) )
  Group = rep( c("Case","Control","Adult","Fetal","Mouse") , c(2,2,3,3,1) )
  
  #####################################################################################
  ##########      get fraction and from fraction to tag since mouse and human coordination different
  #####################################################################################
  
  humanGroup = list( Case=c(1:2), Control=c(3:4), Adult=c(5:7), Fetal=c(8:10)  )
  humanFracs = lapply( humanGroup , function(ind)  getFracFromGffList(humanGffs[ind],humanAnnotation)  )
  humanTags = lapply( humanFracs , function(x)  fracToTag(x) )
  mouseFracs = getFrac(mouseGffs,mouseAnnotation)
  mouseTags = fracToTag(mouseFracs)
  
  Fracs = c( humanFracs , Mouse=list(mouseFracs) )
  Tags = c( humanTags , Mouse=list(mouseTags) )
  
  
  #####################################################################################
  ##########      Venn diagram for all with fetal
  #####################################################################################
  
  
  tagList = Tags[2:5]
  tagAll = unique(do.call(c, tagList))
  
  humanmouseExps0 = lapply( 3:11, function(i) {
    if( Group[i]=="Mouse")
    { 
      fracs = getFrac(Gffs[[i]],mouseAnnotation)
    } else {
      fracs = getFrac(Gffs[[i]],humanAnnotation)
    }
    
    tag = fracToTag(fracs)
    ee = Exps[[i]][ match( tagAll , tag ) ]
    ee[is.na(ee)] = 0
    ee
  } )
  
  humanmouseExps0 = do.call(cbind,humanmouseExps0)
  humanmouseExps = apply(humanmouseExps0,2,function(x) 6000*x/sum(x) )
  exps = t(apply(humanmouseExps,1,function(x) tapply(x,Group[3:11],mean) ) )
  range = which(rowSums(exps[,1:3])>0 & exps[,4]>0)
  hmExp = cbind( log10( rowSums(exps[range,1:3]) + 1 ) , log10(exps[range,4]+1) )
  length(range)
  
  
  #####################################################################################
  ##########      Venn diagram for all with fetal
  #####################################################################################
  
  
  library(VennDiagram)
  indexList <- lapply(tagList, function(x) which(tagAll %in% x))
  names(indexList) = names(tagList)
  setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/comparisonResult/adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100/RevisionThresholdTesting")
  filename = paste0("vennDiagram_withMouse_threeFetals_threshold_",threshold,".png")
  venn.diagram(indexList, filename = filename, height = 2000, 
               width = 2000, fill = 1:length(tagList), alpha = 0.4, 
               imagetype = "png")
  
  
  r = signif( cor.test(hmExp[,2],hmExp[,1])$estimate , 3) 
  p = signif( cor.test(hmExp[,2],hmExp[,1])$p.val , 3) 
  pdf( paste0("Scatterplot_exp_mouse_human_threeFetals_threshold_",threshold,".pdf") )
  plot(hmExp,col='red',main=paste0("r=",r,"; p=",p),xlab='human -log10 read count',ylab='mouse -log10 read count',pch=20,lwd=5)
  abline( lm(hmExp[,2]~hmExp[,1]),lwd=2)
  dev.off()
  

}









