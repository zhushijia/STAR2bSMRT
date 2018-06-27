junc2tag = function(gff) 
{
  lapply( gff , function(read) 
  {
    tmp = paste( read$end[-nrow(read)]+1 , read$start[-1]-1 , sep="," )
    paste(tmp,collapse="; ")
  })
}

s = 50147488
e = 51259000
setwd("/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/gencode/GRCh37.liftover.from.GRCh38.p10.Release26/annotation")
gtf = read.table("NRXN1.gtf",sep="\t")
exon = subset( gtf , V3=='exon' & V4>s & V5<e )
colnames(exon)[c(1,4,5)] = c('chr','start','end')
exon = exon[order(exon$start),]
transcript = sapply( strsplit(as.character(exon$V9),"; ") , function(x) gsub("transcript_id ","",x[2]) )
exon = split( exon[,c(1,4,5)] , transcript )
exon = exon[sapply(exon,nrow)>1]

annotTag = junc2tag(exon)

sum( tag %in% annotTag )
sum( !tag %in% annotTag )



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
  
  exp = lapply( gffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
  }   )
  
  translated = lapply( gffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
  }   )
  
  
  translatedGffs = lapply( 1:length(gffs) , function(i) gffs[[i]][ translated[[i]] ]  )
  names(translatedGffs) = sampleNames
  
  hiPSC = unionGff(translatedGffs[[1]],translatedGffs[[2]])
  humanBrain = unionGff(unionGff(translatedGffs[[6]],translatedGffs[[7]]),translatedGffs[[8]])
  fetal = translatedGffs[[9]]
  patient = unionGff(translatedGffs[[3]],translatedGffs[[4]])
  mergeTranslatedGffs = list( hiPSC=hiPSC , adult=humanBrain , fetal=fetal , patient=patient )
  
  mergeTranslatedTags = lapply( mergeTranslatedGffs , junc2tag )
  
  lapply( mergeTranslatedTags , function(tag) c( sum( tag %in% annotTag ), sum( !tag %in% annotTag )   )   )
  
  
  
}


