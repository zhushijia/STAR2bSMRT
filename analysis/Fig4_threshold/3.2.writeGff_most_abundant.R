library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")
library(Biostrings)
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

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
#samples= c("581","641","2607","553","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks") #
#sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") #
samples= c("581","641","2607","553") #
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2") #
case_ind = 1:2
cont_ind = 3:4
human_annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
case_annotation = human_annotation
cont_annotation = human_annotation
threshold=7

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

fas = lapply( samples , function(sample)
{
    files = dir( paste0(folder,"/",sample,"/",parai) )
    fa = paste0(folder,"/",sample,"/",parai,"/",files[grepl("fa$",files)])
    readDNAStringSet( fa )
} )

translatedFas = lapply( 1:length(fas) , function(i) {
    translatedIds =  sapply( strsplit(names(translatedGffs[[i]])," ") ,function(x) gsub(";", "", x[4])  )
    index = which( names(fas[[i]]) %in% translatedIds )
    fas[[i]][index]  
})

caseGff = NULL
for (i in case_ind) caseGff = unionGff(caseGff,translatedGffs[[i]])
contGff = NULL
for (i in cont_ind) contGff = unionGff(contGff,translatedGffs[[i]])

case_fracs = getFrac(caseGff,case_annotation)
cont_fracs = getFrac(contGff,cont_annotation)

case_fracs2 = as.data.frame(round(case_fracs))
cont_fracs2 = as.data.frame(round(cont_fracs))

range = matchGff(contGff,caseGff)
range = range[!is.na(range)]

caseUniFrac = -case_fracs2[-range, ] # -1 means case unique 
commonFrac = case_fracs2[range, ]

caseUniGff = caseGff[-range]
commonGff = caseGff[range]

range = matchGff(caseGff,contGff)
range = range[!is.na(range)]

contUniFrac = 2*cont_fracs2[-range, ]
contUniGff = contGff[-range]

unique(do.call(c,caseUniFrac))
unique(do.call(c,contUniFrac))
unique(do.call(c,commonFrac))

fracs = rbind(caseUniFrac,contUniFrac,commonFrac)

allGff = c( caseUniGff, contUniGff , commonGff )
exps = sapply( 1:length(translatedGffs) , function(i) {
    ee = translatedExps[[i]][ matchGff( allGff , translatedGffs[[i]] ) ]
    ee[is.na(ee)] = 0
    ee
} )

GffName = sapply( strsplit(names(allGff)," ") ,function(x) gsub(";", "", x[4]) ) 
allFaTmp = do.call(c,translatedFas)
allFa = allFaTmp[ match(GffName,names(allFaTmp)) ]
GffName == names(allFa)

names(allGff) = names(allFa) = apply(exps,1,function(x) paste(paste0(sampleNames,'_exp',x),collapse="_")  )
fracs_ordered = fracs[order( rowSums(exps) , decreasing=T ) ,]
allGff_ordered = allGff[ order( rowSums(exps) , decreasing=T ) ]
exps_ordered = exps[ order( rowSums(exps) , decreasing=T ), ]
allFa_ordered = allFa[ order( rowSums(exps) , decreasing=T ) ]

   
mutantGff = allGff_ordered[rowSums(exps_ordered[,cont_ind]) == 0]
mutantFa = allFa_ordered[rowSums(exps_ordered[,cont_ind]) == 0]
wildtypeGff = allGff_ordered[rowSums(exps_ordered[,cont_ind]) > 0]
wildtypeFa = allFa_ordered[rowSums(exps_ordered[,cont_ind]) > 0]

outputDir = paste0(folder,"/comparisonResult","/",parai,"/gff_fa_thresholded")
setwd(outputDir)
writeXStringSet( mutantFa , "mutantNRXN1.fa" )
writeXStringSet( wildtypeFa , "wildtypeNRXN1.fa" )
writeGff( mutantGff , file = "mutantNRXN1.gff" )
writeGff( wildtypeGff , file = "wildtypeNRXN1.gff" )

