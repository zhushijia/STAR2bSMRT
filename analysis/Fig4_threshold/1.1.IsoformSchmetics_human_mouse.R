gffInfo = function(gff, annotation, threshold)
{
  translated = sapply( strsplit(names(gff),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
  translatedGff = gff[ translated ] 
  translatedExp = sapply( strsplit(names(translatedGff),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
  
  translatedGff = translatedGff[ translatedExp>=threshold ]
  translatedExp = translatedExp[ translatedExp>=threshold ]
  
  translatedFrac = getFrac(translatedGff, annotation)
  translatedFrac2 = as.data.frame(round(translatedFrac))
  translatedTag = apply(translatedFrac2,1,function(x) paste(colnames(translatedFrac2)[x>0],collapse="_") )
  list(gff=translatedGff, exp=translatedExp, tag=translatedTag,frac=translatedFrac2)
}
  
  
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


heatmapMannualColor = function(X,case_col,shared_col)
{
  library(ggplot2)
  library(reshape)
  library(plyr)
  X.m <- melt(X)
  X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
  gg <- ggplot(X.s, aes(x=variable, y=Name))
  gg <- gg + geom_tile(aes(fill = factor(value)), colour = "white")
  #gg <- gg + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred")
  gg <- gg + scale_fill_manual(values=c(case_col,"white",shared_col))
  gg <- gg + labs(x="", y="")
  gg <- gg + theme_bw()
  gg <- gg + theme(panel.grid=element_blank(), panel.border=element_blank())
  base_size <- 9
  gg <- gg + theme(axis.ticks=element_blank(), 
                   axis.text.x=element_text(size=base_size*0.8, angle=300, 
                                            hjust = 0, colour="grey50"))
  print(gg)
  
  
}


heatmapMannualValidateColor = function(X)
{
  library(ggplot2)
  library(reshape)
  library(plyr)
  X.m <- melt(X)
  X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
  gg <- ggplot(X.s, aes(x=variable, y=Name))
  gg <- gg + geom_tile(aes(fill = factor(value)), colour = "lightgrey")
  #gg <- gg + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred")
  gg <- gg + scale_fill_manual(values=c("white","black"))
  gg <- gg + labs(x="", y="")
  gg <- gg + theme_bw()
  gg <- gg + theme(panel.grid=element_blank(), panel.border=element_blank())
  base_size <- 9
  gg <- gg + theme(axis.ticks=element_blank(), 
                   axis.text.x=element_text(size=base_size*0.8, angle=300, 
                                            hjust = 0, colour="grey50"))
  print(gg)
  
  
}


heatmapGradientColor = function(X)
{
  library(ggplot2)
  library(reshape)
  library(plyr)
  X.m <- melt(X)
  X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
  gg <- ggplot(X.s, aes(x=variable, y=Name))
  gg <- gg + geom_tile(aes(fill = value), colour = "white")
  gg <- gg + scale_fill_gradient2(low = "white", high = "darkblue")
  #gg <- gg + scale_fill_manual(values=c("red","white","orange","green"))
  gg <- gg + labs(x="", y="")
  gg <- gg + theme_bw()
  gg <- gg + theme(panel.grid=element_blank(), panel.border=element_blank())
  base_size <- 9
  gg <- gg + theme(axis.ticks=element_blank(), 
                   axis.text.x=element_text(size=base_size*0.8, angle=300, 
                                            hjust = 0, colour="grey50"))
  print(gg)
  
  
}





library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
human_annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
mouse_annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT/NRXN1_mm10_ExonAnnotations_shijia.txt',sep='\t',header=T)

threshold=7

samples= c("581","641","2607","553","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks") #
sampleNames = c("3Del1","3Del2","Cont1","Cont2","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") #

humanGffs = lapply( samples , function(sample)
{
  files = dir( paste0(folder,"/",sample,"/",parai) )
  gff = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
  readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(humanGffs) = sampleNames

humanGffInfo = lapply( humanGffs, function(x) gffInfo( x, human_annotation, 7 )  )
humanTag = unique( do.call(c,lapply(humanGffInfo,function(x)x$tag)) )


folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT"
gff = paste0(folder,"/Exp_SRR1184043_STARlongNew/isoform_ts70_td13.gff" )
mouseGff = readGff( gff , chrom="chr17" , s= 90036900 , e = 91089605 )
mouseGffInfo = gffInfo( mouseGff, mouse_annotation, 7 )
mouseTag = mouseGffInfo$tag
mouseFrac = mouseGffInfo$frac
mouseExp = mouseGffInfo$exp

range = which(mouseTag%in%humanTag)
mouseFrac[-range,] = -mouseFrac[-range,]
Name = reorder( paste0('isoform',1:length(mouseExp)) , mouseExp )
fracs = data.frame( Name , mouseFrac ) 

##########################################################################################
#################   mouse isoform which overlap with human
##########################################################################################

setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT/Exp_SRR1184043_STARlongNew")
pdf("IsoformSchematics_mouseExp.barplot_overlapWithHuman.pdf")
heatmapMannualColor(fracs, case_col="purple", shared_col="darkgreen" )
cols = rep('purple',length(mouseExp) )
cols[range] = 'darkgreen'
barplot(mouseExp[order(Name)],horiz=T,col=cols[order(Name)],border='white')
dev.off()

##########################################################################################
#################   human isoform which overlap with mouse
##########################################################################################

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples= c("581","641","2607","553","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks") #
sampleNames = c("3Del1","3Del2","Cont1","Cont2","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") #

humanGffs = lapply( samples[-c(1:2)] , function(sample)
{
  files = dir( paste0(folder,"/",sample,"/",parai) )
  gff = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
  readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(humanGffs) = sampleNames[-c(1:2)]

humanGffInfo = lapply( humanGffs, function(x) gffInfo( x, human_annotation, 7 )  )
humanTag = unique( do.call(c,lapply(humanGffInfo,function(x)x$tag)) )

humanFrac = do.call( rbind, lapply( humanTag , function(tag) {
  tmp = rep(0,nrow(human_annotation))
  tmp[ as.character(human_annotation$Exon) %in% strsplit(tag,"_")[[1]]  ] = 1
  tmp
} ) )
colnames(humanFrac) = as.character(human_annotation$Exon)

humanExps = do.call( cbind, lapply( humanGffInfo , function(info) {
  tmp = rep(0,length(humanTag))
  tmp[ match( info$tag, humanTag)  ] = info$exp
  tmp
} ) )
normalizedFactor = colSums(humanExps)/min(colSums(humanExps))
humanExps = sapply( 1:ncol(humanExps) , function(i) humanExps[,i]/normalizedFactor[i] )
colnames(humanExps) = names(humanGffs)

range = which(humanTag%in%mouseTag)
humanFrac[-range,] = -humanFrac[-range,]
Name = reorder( paste0('isoform',1:nrow(humanExps)) , rowSums(log10(humanExps+1)) )
fracs = data.frame( Name , humanFrac ) 
logExps = data.frame( Name, log10(humanExps+1) )

setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT/Exp_SRR1184043_STARlongNew")
pdf("IsoformSchematics_humanExp.barplot_overlapWithMouse.pdf")
heatmapMannualColor(fracs, case_col="purple", shared_col="darkgreen" )
heatmapGradientColor(logExps)
cols = rep('purple',nrow(humanExps) )
cols[range] = 'darkgreen'
barplot( rowSums(log10(humanExps+1))[order(Name)],horiz=T,col=cols[order(Name)],border='white')
dev.off()

