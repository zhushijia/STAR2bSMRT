library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")
library(Biostrings)
writeTranslatedGff = function(samples, sampleNames, threshold, outputDir, fileName)
{
  ################################################################################
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
  
  ################################################################################
  
  fas = lapply( samples , function(sample)
  {
    files = dir( paste0(folder,"/",sample,"/",parai) )
    fa = paste0(folder,"/",sample,"/",parai,"/",files[grepl("fa$",files)])
    readDNAStringSet( fa )
  } )
  names(fas) = sampleNames
  
  translatedFas = lapply( 1:length(fas) , function(i) {
    translatedIds =  sapply( strsplit(names(translatedGffs[[i]])," ") ,function(x) gsub(";", "", x[4])  )
    index = which( names(fas[[i]]) %in% translatedIds )
    fas[[i]][index]  
  })
  names(translatedFas) = names(gffs)
  ################################################################################
  setwd(outputDir)
  
  for( i in 1:length(translatedFas) )
  {
    fastaName = paste0( fileName, threshold, "_", names(translatedFas)[i],".fa")
    writeXStringSet( translatedFas[[i]] , fastaName )
    
    gffName = paste0( fileName, threshold, "_", names(translatedFas)[i],".gff")
    writeGff( isoform=translatedGffs[[i]] , file = gffName )
  }
  
}

##################################
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples= c("581","641","2607","553","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks") #
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") #
threshold=7
fileName = "TranslatedThresholed"
outputDir = paste0(folder,"/comparisonResult","/",parai,"/gff_fa_thresholded")
writeTranslatedGff(samples, sampleNames, threshold, outputDir, fileName)
  
##################################
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples = c("553_GABA","553_NGN2","553")
sampleNames = c("553_GABA","553_NGN2","553_FBN")
threshold=7
fileName = "TranslatedThresholed"
outputDir = paste0(folder,"/comparisonResult","/",parai,"/gff_fa_thresholded")
writeTranslatedGff(samples, sampleNames, threshold, outputDir, fileName)

