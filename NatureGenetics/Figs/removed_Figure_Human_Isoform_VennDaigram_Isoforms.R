
drawGffVennDiagram = function(samples, threshold, groups, outputDir, fileName)
{
  human_nrxn1_chr = "chr2"
  human_nrxn1_start = 50149082 
  human_nrxn1_end = 51255411
  
  gffs = lapply( samples , function(sample)
  {
    gff = paste0("Data/Human_NRXN1alpha_Translated_Thres7_",sample,".gff")
    readGff( gff , chrom=human_nrxn1_chr , s=human_nrxn1_start , e=human_nrxn1_end )
  } )
  names(gffs) = samples

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

##########################################################################################
#################   load codes for heatmap and gff
##########################################################################################

source("SourceCode_Gff_Op.R")
source("SourceCode_Heatmap.R")

##########################################################################################
#################   human samples and basic information
##########################################################################################
ids = c("581","641","2607","553","642",
        "NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13",
        "fetal","fetal_23wks","fetal_3wks",
        "553_GABA","553_NGN2","553") #

samples = c("p3Del1","p3Del2","Cont1","Cont2","Cont3",
            "Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3",
            "553_GABA","553_NGN2","553_FBN") #
threshold = 7
outputDir = getwd()
fileName = "VennDiagram"

##########################################################################################
#################   venn diagram for case, cont, adult and fetal
##########################################################################################

groups = list( Case=c(1,2), Cont=c(3:5), Adult=c(6:8), Fetal=c(9:11) ) 
drawGffVennDiagram(samples, threshold, groups, outputDir, fileName)
  
##########################################################################################
#################   venn diagram for cont, adult and fetal
##########################################################################################

groups = list( Cont=c(3:5), Adult=c(6:8), Fetal=c(9:11) ) 
drawGffVennDiagram(samples, threshold, groups, outputDir, fileName)

##########################################################################################
#################   venn diagram for case, and cont
##########################################################################################

groups = list( Case=c(1,2), Cont=c(3:4) ) 
drawGffVennDiagram(samples, threshold, groups, outputDir, fileName)

##########################################################################################
#################   venn diagram for gaba, ngn2, and forebrain neuron
##########################################################################################

groups = list( GABA=c(12), NGN2=c(13), FBN=c(14) ) 
drawGffVennDiagram(samples, threshold, groups, outputDir, fileName)

