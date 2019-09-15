##########################################################################################
#################   load codes for heatmap and gff
##########################################################################################

source("SourceCode_Gff_Op.R")
source("SourceCode_Heatmap.R")

##########################################################################################
#################   load basic information
##########################################################################################

human_annotation = read.table('Data/human_NRXN1_Exon_Annotation_hg19.txt',sep='\t',header=T)
mouse_annotation = read.table('Data/mouse_NRXN1_Exon_Annotation_mm10.txt',sep='\t',header=T)

human_nrxn1_chr = "chr2"
human_nrxn1_start = 50149082 
human_nrxn1_end = 51255411

mouse_nrxn1_chr = "chr17"
mouse_nrxn1_start = 90036900
mouse_nrxn1_end = 91089605

##########################################################################################
#################   human samples 
##########################################################################################

samples= c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") 

ids = c("581","641","2607","553","642",
        "NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13",
        "fetal","fetal_23wks","fetal_3wks",
        "553_GABA","553_NGN2","553") #

samples = c("p3Del1","p3Del2","Cont1","Cont2","Cont3",
            "Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3",
            "553_GABA","553_NGN2","553_FBN") #

humanGffs = lapply( samples , function(sample)
{
  gff = paste0("Data/Human_NRXN1alpha_Translated_Thres7_",sample,".gff")
  readGff( gff , chrom=human_nrxn1_chr , s=human_nrxn1_start , e=human_nrxn1_end )
} )
names(humanGffs) = samples
sapply(humanGffs,length)

##########################################################################################
#################   mouse samples 
##########################################################################################

mouseGff = readGff( "Data/mouse_ts70_td13.gff" , chrom=mouse_nrxn1_chr , s=mouse_nrxn1_start , e=mouse_nrxn1_end )
mouseGffInfo = gffInfo( gff=mouseGff, annotation=mouse_annotation, threshold=7 )
mouseTag = mouseGffInfo$tag


library(VennDiagram)
dir.create("Figs")

##########################################################################################
#################   draw venn diagram
##########################################################################################

groups = list( Adult=c(6:8), Fetal=c(9:11) ) 
groupGffs = lapply(groups, function(group) {
    tmp = NULL
    for (i in group) tmp = unionGff(tmp,humanGffs[[i]])
    tmp
})
names(groupGffs) = names(groups)
groupGffInfo = lapply( groupGffs, function(x) gffInfo( gff=x, annotation=human_annotation, threshold=7 )  )
groupTag = lapply(groupGffInfo,function(x)x$tag)
groupTag$mouse = mouseTag

tagAll = unique(do.call(c, groupTag))
indexList <- lapply(groupTag, function(x) which(tagAll %in% x))
size = sapply(indexList,length)
names(indexList) = paste0( names(groupTag), "(", size , ")")
filename = "Figs/Fig2A_vennDiagram_Mouse_and_Human.png"
venn.diagram(indexList, filename = filename, height = 2000, 
             width = 2000, fill = 1:length(groupTag), alpha = 0.4, 
             imagetype = "png")


##########################################################################################
#################   draw venn diagram
##########################################################################################

groups = list( Cont=c(3:4), Adult=c(6:8), Fetal=c(9:11) ) 
groupGffs = lapply(groups, function(group) {
    tmp = NULL
    for (i in group) tmp = unionGff(tmp,humanGffs[[i]])
    tmp
})
names(groupGffs) = names(groups)
groupGffInfo = lapply( groupGffs, function(x) gffInfo( gff=x, annotation=human_annotation, threshold=7 )  )
groupTag = lapply(groupGffInfo,function(x)x$tag)

tagAll = unique(do.call(c, groupTag))
indexList <- lapply(groupTag, function(x) which(tagAll %in% x))
size = sapply(indexList,length)
names(indexList) = paste0( names(groupTag), "(", size , ")")
filename = "Figs/Fig2E_vennDiagram_Human_hiPSC_Adult_Fetal.png"
venn.diagram(indexList, filename = filename, height = 2000, 
             width = 2000, fill = 1:length(groupTag), alpha = 0.4, 
             imagetype = "png")

##########################################################################################
#################   draw venn diagram
##########################################################################################
groups = list( GABA=c(12), NGN2=c(13), FBN=c(14) ) 
groupGffs = lapply(groups, function(group) {
    tmp = NULL
    for (i in group) tmp = unionGff(tmp,humanGffs[[i]])
    tmp
})
names(groupGffs) = names(groups)
groupGffInfo = lapply( groupGffs, function(x) gffInfo( gff=x, annotation=human_annotation, threshold=7 )  )
groupTag = lapply(groupGffInfo,function(x)x$tag)

tagAll = unique(do.call(c, groupTag))
indexList <- lapply(groupTag, function(x) which(tagAll %in% x))
size = sapply(indexList,length)
names(indexList) = paste0( names(groupTag), "(", size , ")")
filename = "Figs/Fig3A_vennDiagram_Human_GABA_NGN2_FBN.png"
venn.diagram(indexList, filename = filename, height = 2000, 
             width = 2000, fill = 1:length(groupTag), alpha = 0.4, 
             imagetype = "png")

