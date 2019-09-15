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

samples= c("p3Del1","p3Del2","Cont1","Cont2","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") 

humanGffs = lapply( samples , function(sample)
{
  gff = paste0("Data/Human_NRXN1alpha_Translated_Thres7_",sample,".gff")
  readGff( gff , chrom=human_nrxn1_chr , s=human_nrxn1_start , e=human_nrxn1_end )
} )
names(humanGffs) = samples
sapply(humanGffs,length)

humanGffInfo = lapply( humanGffs, function(x) gffInfo( gff=x, annotation=human_annotation, threshold=7 )  )
humanTag = unique( do.call(c,lapply(humanGffInfo,function(x)x$tag)) )

##########################################################################################
#################   mouse samples 
##########################################################################################

mouseGff = readGff( "Data/mouse_ts70_td13.gff" , chrom=mouse_nrxn1_chr , s=mouse_nrxn1_start , e=mouse_nrxn1_end )
mouseGffInfo = gffInfo( gff=mouseGff, annotation=mouse_annotation, threshold=7 )
mouseTag = mouseGffInfo$tag
mouseFrac = mouseGffInfo$frac
mouseExp = mouseGffInfo$exp

range = which(mouseTag%in%humanTag)
mouseFrac[-range,] = -mouseFrac[-range,]
Name = reorder( paste0('isoform',1:length(mouseExp)) , mouseExp )
fracs = data.frame( Name , mouseFrac ) 

##########################################################################################
#################   draw schematics
##########################################################################################
dir.create("Figs")

pdf("Figs/Fig2C_Mouse_Isoform_Schematics_overlapped_with_human.pdf")
# draw isoform schematics
heatmapMannualColor2(fracs, case_col="purple", shared_col="darkgreen" )

# draw expression barplot
cols = rep('purple',length(mouseExp) )
cols[range] = 'darkgreen'
barplot(mouseExp[order(Name)],horiz=T,col=cols[order(Name)],border='white')

dev.off()




