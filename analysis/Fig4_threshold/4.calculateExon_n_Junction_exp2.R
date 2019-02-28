library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")


heatmapGradientColor = function(X,Y,title)
{
  library(ggplot2)
  library(reshape)
  library(plyr)
  library(gridExtra)
  library(grid)
  library(lattice)
  
  X.m <- melt(X)
  X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
  X.s$Name <- ordered( X.s$Name , rev( levels(X.s$Name) ) )
  gg <- ggplot(X.s, aes(x=variable, y=Name))
  gg <- gg + geom_tile(aes(fill = value), colour = "white")
  gg <- gg + scale_fill_gradient2(low = "white", high = "darkblue")
  #gg <- gg + scale_fill_manual(values=c("red","white","orange","green"))
  gg <- gg + labs(x="", y="")
  gg <- gg + theme_bw()
  gg <- gg + theme(panel.border=element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks.x = element_blank()) #
  base_size <- 9
  gg <- gg + theme(axis.ticks=element_blank(), 
                   axis.text.x=element_text(size=base_size*0.8, angle=300, 
                                            hjust = 0, colour="grey50"),
                   legend.position = "bottom", 
                   legend.direction = "horizontal")
  
  p <- ggplot(data=Y, aes(x=Name, y=y)) + geom_bar(stat="identity") +
    ggtitle(title)
  p.nox <- p + theme(legend.position = "none",
                          axis.title = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks.x = element_blank())
  
  grid.arrange(p.nox, gg, nrow=2, heights=0.7:2.5)
  
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

calculateExonExp = function(cc_fracs,ccExp)
{
  exonExpExpEachTranscript = apply( cc_fracs , 2 , function(x) x*rowSums(ccExp) )
  colSums(exonExpExpEachTranscript)
}

calculateJuncFrac = function(cc_fracs,ccExp)
{
  isoform = apply(cc_fracs , 1 , function(x) {
    iso1 = paste( which(x>0), collapse = "_" )
    paste0( "_" , iso1 , "_")
  })
  
  
  junc = matrix(ncol=ncol(cc_fracs),nrow=ncol(cc_fracs),data=0)
  colnames(junc) = colnames(cc_fracs)
  rownames(junc) = colnames(cc_fracs)
  
  for(i in 1:ncol(cc_fracs))
  {
    for(j in i:ncol(cc_fracs))
    {
      if(i==j)
      {
        junc[i,i] = sum( grepl( paste("",i,"",sep="_") , isoform ) )
      } else {
        junc[i,j] = sum( grepl( paste("",i,j,sep="_") , isoform ) )
      }
    }
    
  }
  
  frac = junc/diag(junc)
  frac[lower.tri(frac)] = -1
  diag(frac) = 0
  frac
}



drawJuncExp = function(samples, sampleNames, threshold, groups, outputDir, fileName)
{
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
  
  allGffs = NULL
  for (i in do.call(c,groups) ) allGffs = unionGff(allGffs,translatedGffs[[i]])
  
  groupGffs = lapply(groups, function(group) {
    tmp = NULL
    for (i in group) tmp = unionGff(tmp,translatedGffs[[i]])
    tmp
  })
  names(groupGffs) = names(groups)
  
  groupFracs = lapply( groupGffs, function(gff) {
    frac = getFrac(gff,annotation)
    as.data.frame(round(frac))
  } )
  
  groupExps = lapply( 1:length(groups) , function(gp) {
    sapply( groups[[gp]], function(i) {
      ee = translatedExps[[i]][ matchGff( groupGffs[[gp]] , translatedGffs[[i]] ) ]
      ee[is.na(ee)] = 0
      ee })
  } )
    
  groupExonExp = lapply( 1:length(groups) , function(i) 
    calculateExonExp(groupFracs[[i]],groupExps[[i]]) )
  names(groupExonExp) = names(groups)
  
  groupJuncFrac = lapply( 1:length(groups) , function(i) 
    calculateJuncFrac(groupFracs[[i]],groupExps[[i]]) )
  names(groupJuncFrac) = names(groups)
 
  library(gplots)
  setwd(outputDir)
  pdf( paste0(fileName,"_1.pdf"), w=3,h=5 )
  for(i in 1:length(groups))
  {
    x = groupJuncFrac[[i]]
    x[is.na(x)] = 0
    Name = reorder( rownames(x) , 1:nrow(x) )
    X = data.frame(Name,x)
    Y = data.frame(Name,y=groupExonExp[[i]])
    heatmapGradientColor(X,Y,names(groups)[i])
  }
  dev.off()
  
  pdf( paste0(fileName,"_2.pdf") )
  for(i in 1:length(groups))
  {
    barplot( groupExonExp[[i]], border=F, las=3, main=names(groups)[i] )
    color = list( redgreen(256) , topo.colors(256) , cm.colors(256) , terrain.colors(256) , rainbow(256), heat.colors(256)  )
    heatmap.2(groupJuncFrac[[i]], main=names(groups)[i], col=color[[3]], scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
  }
  dev.off()
  
  
}


annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
outputDir = paste0(folder,"/comparisonResult","/",parai)
samples= c("581","641","2607","553","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks","553_GABA","553_NGN2") #
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3","553_GABA","553_NGN2") #
groups = as.list(1:length(samples))
names(groups) = sampleNames
threshold=7
fileName = "groups_junc_exp"

drawJuncExp(samples, sampleNames, threshold, groups, outputDir, fileName)


annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
outputDir = paste0(folder,"/comparisonResult","/",parai)
samples= c("581","641","2607","553","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks","553_GABA","553_NGN2") #
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3","553_GABA","553_NGN2") #
groups = list( Case=c(1,2), Cont=c(3:5), Fetal=c(6:8), Adult=c(9:11), 12, 13 ) 
names(groups) = c("Case","Control","Fetal","Adult","553_GABA","553_NGN2")
threshold=7
fileName = "groups_case_cont_junc_exp"

drawJuncExp(samples, sampleNames, threshold, groups, outputDir, fileName)





