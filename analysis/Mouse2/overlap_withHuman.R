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




humanAnnotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
mouseAnnotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT/NRXN1_mm10_ExonAnnotations_shijia.txt',sep='\t',header=T)
mouseAnnotation$Start[1] = 91088726
mouseAnnotation$Stop[26] = 90037077

library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")

samples= c("581","641","2607","553","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal")
sampleNames = c("3' Del 1","3' Del 2","Control 1","Control 2","Adult 1","Adult 2","Adult 3","Fetal")

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai = "adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"

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


caseGff = unionGff(translatedGffs[[1]],translatedGffs[[2]])
contGff = unionGff(translatedGffs[[3]],translatedGffs[[4]])
adultGff = unionGff(translatedGffs[[5]],translatedGffs[[6]])
adultGff = unionGff(adultGff,translatedGffs[[7]])
fetalGff = translatedGffs[[8]]

case_fracs = getFrac(caseGff,humanAnnotation)
cont_fracs = getFrac(contGff,humanAnnotation)
adult_fracs = getFrac(adultGff,humanAnnotation)
fetal_fracs = getFrac(fetalGff,humanAnnotation)



gff = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT/Exp_SRR1184043_STARlongNew/isoform_ts70_td13.gff"
mouseGff = readGff( gff , chrom="chr17" , s= 90036900 , e = 91089605 )
translated = sapply( strsplit(names(mouseGff),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
mouseGff = mouseGff[translated]
mouse_exps = sapply( strsplit(names(mouseGff),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
mouse_fracs = getFrac(mouseGff,mouseAnnotation)

caseIsoform = fracToTag(case_fracs)
contIsoform = fracToTag(cont_fracs)
adultIsoform = fracToTag(adult_fracs)
fetalIsoform = fracToTag(fetal_fracs)
mouseIsoform = fracToTag(mouse_fracs)



setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/comparisonResult/adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100")

tagList = list(case=caseIsoform,cont=contIsoform,adult=adultIsoform,fetal=fetalIsoform,mouse=mouseIsoform)
library(VennDiagram)
tagAll = unique(do.call(c, tagList))
indexList <- lapply(tagList, function(x) which(tagAll %in% x))
names(indexList) = names(tagList)
filename = "vennDiagram_withMouse_withCase.png"
venn.diagram(indexList, filename = filename, height = 2000, 
             width = 2000, fill = 1:length(tagList), alpha = 0.4, 
             imagetype = "png")


tagList = list(cont=contIsoform,mouse=mouseIsoform,fetal=fetalIsoform,adult=adultIsoform)
library(VennDiagram)
tagAll = unique(do.call(c, tagList))
indexList <- lapply(tagList, function(x) which(tagAll %in% x))
names(indexList) = names(tagList)
filename = "vennDiagram_withMouse.png"
venn.diagram(indexList, filename = filename, height = 2000, 
             width = 2000, fill = 1:length(tagList), alpha = 0.4, 
             imagetype = "png")

humanTagList = unique( do.call(c, tagList[-2]) )
mo = mouse_exps[ tagList$mouse %in% humanTagList ]
mn = mouse_exps[ !tagList$mouse %in% humanTagList ]

pdf('mouse_exp_overlapWithHuman.pdf')
boxplot( list(mo , mn) , col=c('darkred','darkgreen') , outline=F  )
dev.off()

tagList4 = tagList
tagAll4 = tagAll

#####################################################################################

translatedHumanExps = lapply( 3:8, function(i) {
  gff = translatedGffs[[i]]
  fracs = getFrac(gff,humanAnnotation)
  tag = fracToTag(fracs)
  ee = translatedExps[[i]][ match( tagAll4 , tag ) ]
  ee[is.na(ee)] = 0
  ee
} )

translatedHumanExps = do.call(cbind,translatedHumanExps)
translatedHumanExps = apply(translatedHumanExps,2,function(x) 6000*x/sum(x) )
colnames(translatedHumanExps) = c("Control","Control","Adult","Adult","Adult","Fetal")
HumanExp = apply(translatedHumanExps,1,function(x) sum( tapply(x,colnames(translatedHumanExps),mean)) )

MouseExp = mouse_exps[ match( tagAll4 , tagList$mouse ) ]
MouseExp[is.na(MouseExp)] = 0


hmExp = data.frame(log10(HumanExp+1),log10(MouseExp+1))
hmExp = hmExp[ hmExp[,1]>0 & hmExp[,2]>0,]
pdf('Scatterplot_exp_mouse_human.pdf')
plot(hmExp,col='red',xlab='human -log10 read count',ylab='mouse -log10 read count',pch=20,lwd=5)
abline( lm(hmExp[,2]~hmExp[,1]),lwd=2)
dev.off()


cor.test(hmExp[,2],hmExp[,1])










