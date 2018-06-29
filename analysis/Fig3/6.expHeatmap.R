heatmapMannualColor = function(X)
{
  library(ggplot2)
  library(reshape)
  library(plyr)
  pdf("case.cont.heatmap_translated_ConsiderOverlap.pdf")
  X.m <- melt(X)
  X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
  gg <- ggplot(X.s, aes(x=variable, y=Name))
  gg <- gg + geom_tile(aes(fill = factor(value)), colour = "white")
  #gg <- gg + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred")
  gg <- gg + scale_fill_manual(values=c("red","white","orange","green"))
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

heatmapGradientColor = function(X)
{
  library(ggplot2)
  library(reshape)
  library(plyr)
  pdf("case.cont.heatmap_translated_Exp.pdf")
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
  dev.off()
  
  
}


library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

#annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1.txt',sep='\t')
annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1ExonAnnotations.txt',sep='\t',header=T)
annotation[2,2] = annotation[3,3]+1
annotation[25,3] = annotation[24,2]-1

######################################
annotation[27,] = annotation[8,]
annotation[8,3] = 50848363
annotation[27,2] = 50848364
annotation[,1] = as.character(annotation[,1])
annotation[27,1] = "exon7c"
annotation[27,4] = 23

library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone4/setup")

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
samples= c("581","641","2607","553")#,"642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal")
sampleNames = c("3' Del 1","3' Del 2","Control 1","Control 2")#,"Control 3","Adult 1","Adult 2","Adult 3","Fetal")


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

	gffs = gffs[1:4]
	
	translated = lapply( gffs , function(x) {
	sapply( strsplit(names(x),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
	}   )

	translatedGffs = lapply( 1:length(gffs) , function(i) gffs[[i]][ translated[[i]] ]  )
	names(translatedGffs) = names(gffs)

	translatedExps = lapply( translatedGffs , function(x) {
	  sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
	}   )
	names(translatedExps) = names(gffs)
	
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
  	fracs
  	
	}

	caseGff = unionGff(translatedGffs[[1]],translatedGffs[[2]])
	contGff = unionGff(translatedGffs[[3]],translatedGffs[[4]])
	case_fracs = getFrac(caseGff,annotation)
	cont_fracs = getFrac(contGff,annotation)

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
	
	#allGff = unionGff( unionGff( caseUniGff, contUniGff ) , commonGff )
	allGff = c( caseUniGff, contUniGff , commonGff )
	#allGff = allGff[ match( rownames(fracs),names(allGff)) ]
	exps = sapply( 1:length(translatedGffs) , function(i) {
	  ee = translatedExps[[i]][ matchGff( allGff , translatedGffs[[i]] ) ]
	  ee[is.na(ee)] = 0
	  ee
	} )
  colnames(exps) = names(translatedGffs)
	exps = data.frame(exps)
	
	cased = matchGff(caseUniGff,allGff)
	contd = matchGff(contUniGff,allGff)
	shared = matchGff(commonGff,allGff)
	
	
	cols = rep('',nrow(exps))
	cols[shared] = "orange"
	cols[cased] = "red"
	cols[contd] = "green"
	
	de = t( apply( exps , 1 , function(x) {
	  case = mean(x[1:2])
	  cont = mean(x[3:4])
	  fc = log2(case+1)-log2(cont+1)
	} ) )
	
	de = as.numeric(de)
	
	Name = reorder( paste0('isoform',1:nrow(exps)) , de )
	logExps = data.frame( Name, log10(exps+1) )
	fracs = data.frame( Name , fracs ) 
	
	
	
	outputDir = paste0(folder,"/comparisonResult","/",parai)
	system( paste( "mkdir -p" , outputDir ) )
	setwd(outputDir)



	heatmapMannualColor(fracs)
	heatmapGradientColor(logExps)
	
	pdf("barplot.pdf",h=10,w=4)
	barplot(de[order(Name)],col=cols[order(Name)],horiz=T,border='grey')
	dev.off()
	
	
}
