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


heatmapMannualColor = function(X)
{
  library(ggplot2)
  library(reshape)
  library(plyr)
  pdf("mouse.heatmap_isoforms_overlapWithHuman.pdf")
  X.m <- melt(X)
  X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
  gg <- ggplot(X.s, aes(x=variable, y=Name))
  gg <- gg + geom_tile(aes(fill = factor(value)), colour = "white")
  #gg <- gg + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred")
  gg <- gg + scale_fill_manual(values=c("darkred","white","darkgreen"))
  gg <- gg + labs(x="", y="")
  gg <- gg + theme_bw()
  gg <- gg + theme(panel.grid=element_blank(), panel.border=element_blank())
  base_size <- 9
  gg <- gg + theme(axis.ticks=element_blank(), 
             axis.text.x=element_text(size=base_size*0.8, angle=300, 
                                      hjust = 0, colour="grey50"))
  print(gg)
  dev.off()
  
  
}


# this one is first loaded from the code "overlap_withHuman.R"
range = which( tagList$mouse %in% humanTagList )

annotation = read.table('/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT/NRXN1_mm10_ExonAnnotations_shijia.txt',sep='\t',header=T)

library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")

folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT"
gff = paste0(folder,"/Exp_SRR1184043_STARlongNew/isoform_ts70_td13.gff" )
mouseGff = readGff( gff , chrom="chr17" , s= 90036900 , e = 91089605 )
translated = sapply( strsplit(names(mouseGff),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
mouseGff = mouseGff[translated]
exps = sapply( strsplit(names(mouseGff),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
mouse_fracs = getFrac(mouseGff,annotation)
mouse_fracs[range,] = -mouse_fracs[range,]
Name = reorder( paste0('isoform',1:length(exps)) , exps )
fracs = data.frame( Name , round(mouse_fracs) ) 


setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT/Exp_SRR1184043_STARlongNew")

################################################################################
#############  isoform schematics
################################################################################

heatmapMannualColor(fracs)

cols = rep('darkgreen',length(exps) )
cols[range] = 'darkred'
pdf("mouseExp.barplot_overlapWithHuman.pdf",h=10,w=4)
barplot(exps[order(Name)],horiz=T,col=cols[order(Name)],border='white')
dev.off()

