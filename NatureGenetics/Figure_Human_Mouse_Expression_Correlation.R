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

ids = c("2607","553",
        "NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13",
        "fetal","fetal_23wks","fetal_3wks") #

samples = c("Cont1","Cont2","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3","Mouse")

Info = lapply( samples , function(sample)
{
    if(sample=='Mouse')
    {
        gff = readGff( "Data/mouse_ts70_td13.gff" , chrom=mouse_nrxn1_chr , s=mouse_nrxn1_start , e=mouse_nrxn1_end )
        info = gffInfo( gff=gff, annotation=mouse_annotation, threshold=7 )
    } else {
        fileName = paste0("Data/Human_NRXN1alpha_Translated_Thres7_",sample,".gff")
        gff = readGff( fileName , chrom=human_nrxn1_chr , s=human_nrxn1_start , e=human_nrxn1_end )
        info = gffInfo( gff=gff, annotation=human_annotation, threshold=7 )
    }
    info
} )
Tags = lapply(Info, function(info) info$tag)
Exps = lapply(Info, function(info) info$exp)


tagAll = unique(do.call(c, Tags))

Exps_mat = lapply( 1:length(Tags), function(i) {
    tag = Tags[[i]]
    ee = Exps[[i]][ match( tagAll , tag ) ]
    ee[is.na(ee)] = 0
    ee
} )

Exps_mat = do.call(cbind,Exps_mat)
Exps_mat = apply(Exps_mat,2,function(x) 6000*x/sum(x) )
Group = rep( c("Control","Adult","Fetal","Mouse") , c(2,3,3,1) )
exps = t(apply(Exps_mat,1,function(x) tapply(x,Group,mean) ) )
range = which(rowSums(exps[,1:3])>0 & exps[,4]>0)
hmExp = cbind( log10( rowSums(exps[range,1:3]) + 1 ) , log10(exps[range,4]+1) )
r = signif( cor.test(hmExp[,2],hmExp[,1])$estimate , 3) 
p = signif( cor.test(hmExp[,2],hmExp[,1])$p.val , 3) 

pdf("Figs/Fig2B_Human_Mouse_Expression_Correlation.pdf",h=5.7,w=5.6)
plot(hmExp,col='darkgreen',main=paste0("r=",r,"; p=",p),xlab='Human log10Read Count',ylab='Mouse log10Read Count',pch=20,lwd=5)
abline( lm(hmExp[,2]~hmExp[,1]),lwd=2)
dev.off()





