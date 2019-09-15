
writeTranslatedGff = function(samples, sampleNames, threshold, outputDir, fileName)
{
  ################################################################################
  gffs = lapply( samples , function(sample)
  {
      loc = paste0(folder,"/",sample,"/STAR2bSMRT")
      gff = dir(loc)[grep("gff",dir(loc))]
      gff = paste0(loc,"/",gff)
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
  
  for( i in 1:length(translatedGffs) ) 
  {
      names(translatedGffs[[i]]) = gsub("gene_id SS.1; transcript_id |;", "", names(translatedGffs[[i]]) )
  }
  ################################################################################
  
  fas = lapply( samples , function(sample)
  {
      loc = paste0(folder,"/",sample,"/STAR2bSMRT")
      fa = dir(loc)[grep("fa",dir(loc))]
      fa = paste0(loc,"/",fa)
      readDNAStringSet( fa )
  } )
  names(fas) = sampleNames
  
  translatedFas = lapply( 1:length(fas) , function(i) {
    index = which( names(fas[[i]]) %in% names(translatedGffs[[i]]) )
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

library(STAR2bSMRT,lib.loc="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/forRelease/r3")
library(Biostrings)
folder="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/test/PackageReleaseResult5"

samples= c("Case1_581","Case2_641","Cont1_2607","Cont2_553","Cont3_642",
           "Adult1_NRXN_adult_dlPFC1_12","Adult2_adult_dlPFC1_10","Adult3_adult_dlPFC1_13",
           "Fetal1_fetal","Fetal2_fetal_23wks","Fetal3_fetal_3wks",
           "GABA_553","NGN2_553","Cont2_553") #
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3","553_GABA","553_NGN2","553_FBN") #
threshold=7
fileName = "Human_NRXN1alpha_Translated_Thres"
outputDir = paste0(folder,"/STAR2bSMRT_gff_fa_thresholded")
writeTranslatedGff(samples, sampleNames, threshold, outputDir, fileName)

