library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")

drawGffVennDiagram = function(samples, sampleNames, threshold, groups, outputDir, fileName)
{
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
  names(translatedGffs) = names(gffs)
  
  translatedExps = lapply( translatedGffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
  }   )
  names(translatedExps) = names(gffs)
  
  translatedGffs = lapply( 1:length(translatedGffs) , function(i) translatedGffs[[i]][ translatedExps[[i]]>=threshold ] )
  translatedExps = lapply( 1:length(translatedGffs) , function(i) translatedExps[[i]][ translatedExps[[i]]>=threshold ] )
  names(translatedGffs) = names(gffs)
  names(translatedExps) = names(gffs)
  
  allGffs = NULL
  for (i in do.call(c,groups) ) allGffs = unionGff(allGffs,translatedGffs[[i]])
  
  groupGffs = lapply(groups, function(group) {
    tmp = NULL
    for (i in group) tmp = unionGff(tmp,translatedGffs[[i]])
    tmp
  })
  names(groupGffs) = names(groups)
  
  setwd(outputDir)
  
  pngName = paste( c(fileName, 
                     paste( names(groups),sapply(groupGffs,length) , sep="" ) ,
                     paste( "ALL", length(allGffs), sep="" ),
                     "threshold", threshold,".png"), collapse="_" )
  
  vennDiagramGff( groupGffs , filename=pngName )
  
}

##################################
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples= c("581","641","2607","553","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks","553_GABA","553_NGN2") #
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3","553_GABA","553_NGN2") #
groups = list( case=c(1:2), AllHumanNoCase=c(3:13) ) 
threshold=7
fileName = "VennDiagram"
outputDir = paste0(folder,"/comparisonResult","/",parai)
drawGffVennDiagram(samples, sampleNames, threshold, groups, outputDir, fileName)


##################################
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples= c("581","641","2607","553","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks") #
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") #
groups = list( Case=c(1,2), Cont=c(3:5), Adult=c(6:8), Fetal=c(9:11) ) 
threshold=7
fileName = "VennDiagram"
outputDir = paste0(folder,"/comparisonResult","/",parai)
drawGffVennDiagram(samples, sampleNames, threshold, groups, outputDir, fileName)
  
##################################
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples= c("581","641","2607","553","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks") #
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") #
groups = list( Cont=c(3:5), Adult=c(6:8), Fetal=c(9:11) ) 
threshold=7
fileName = "VennDiagram"
outputDir = paste0(folder,"/comparisonResult","/",parai)
drawGffVennDiagram(samples, sampleNames, threshold, groups, outputDir, fileName)

##################################
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples= c("581","641","2607","553","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks") #
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") #
groups = list( Case=c(1,2), Cont=c(3:4) ) 
threshold=7
fileName = "VennDiagram"
outputDir = paste0(folder,"/comparisonResult","/",parai)
drawGffVennDiagram(samples, sampleNames, threshold, groups, outputDir, fileName)

##################################
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples = c("553_GABA","553_NGN2","553")
sampleNames = c("553_GABA","553_NGN2","553_FBN")
groups = list( GABA=c(1), NGN2=c(2), FBN=c(3) ) 
threshold=7
fileName = "VennDiagram"
outputDir = paste0(folder,"/comparisonResult","/",parai)
drawGffVennDiagram(samples, sampleNames, threshold, groups, outputDir, fileName)

