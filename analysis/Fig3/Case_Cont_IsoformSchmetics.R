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

heatmapMannualValidateColor = function(X)
{
  library(ggplot2)
  library(reshape)
  library(plyr)
  pdf("case.cont.heatmap_validated_byAdult.pdf")
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





plotExonCoexpress = function(case_fracs, cont_fracs )
{
  
  f = function(X)
  {
    X.m <- melt(X)
    X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
    gg <- ggplot(X.s, aes(x=variable, y=Name))
    gg <- gg + geom_tile(aes(fill = value), colour = "white")
    gg <- gg + scale_fill_gradient2(low = "white", high = "steelblue")
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
 
  library(ggplot2)
  library(reshape)
  library(plyr)
  
  case_dist = as.matrix( dist(t(case_fracs)) )
  cont_dist = as.matrix( dist(t(cont_fracs)) )
  #diff_dist = log(case_dist+1) - log( cont_dist +1 )
  diff_dist = case_dist- cont_dist
  Name = reorder(colnames(case_dist),1:nrow(case_dist))
  case_dist = data.frame(Name,case_dist)
  cont_dist = data.frame(Name,cont_dist)
  diff_dist = data.frame(Name,diff_dist)
  
  pdf("case.cont.Exon_CoExpression.pdf")
  f(case_dist)
  f(cont_dist)
  f(diff_dist)
  dev.off()
  
}

library(DEGseq)	
	
DEGSeq.callDE <- function(  countData , condition , countThres=2 ) 
{
	countData       <- countData[ rowSums(countData) >= countThres , ]  # filtering 
	geneID = sapply(strsplit(rownames(countData),"[.]"),function(x)x[1])
	data = data.frame( geneID=geneID , countData )
	
	uni_condition = unique(condition)
	c1_index = which( condition %in% uni_condition[1] ) + 1
	c2_index = which( condition %in% uni_condition[2] ) + 1
	
	DEGexp(geneExpMatrix1=data, geneCol1=1, expCol1=c1_index, groupLabel1=uni_condition[1], 
		   geneExpMatrix2=data, geneCol2=1, expCol2=c2_index, groupLabel2=uni_condition[1], 
		   method="MARS",outputDir='DEGseq_result')

}		   




annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)

library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")

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
  
	# normalized to sum
	rawExps = exps
	colnames(rawExps) = names(translatedGffs)
	
	normalizedFactor = colSums(exps)/min(colSums(exps))
	exps = sapply( 1:ncol(exps) , function(i) exps[,i]/normalizedFactor[i] )
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

	################################################################################
	#############  isoform schematics
	################################################################################
	
	heatmapMannualColor(fracs)
	heatmapGradientColor(logExps)
	plotExonCoexpress(case_fracs,cont_fracs)
	
	pdf("case.cont.barplot.pdf",h=10,w=4)
	barplot(de[order(Name)],col=cols[order(Name)],horiz=T,border='white')
	dev.off()
	
	pdf("case.cont.pairs.pdf")
	pairs(log10(rawExps+1),pch=16,col='red')
	dev.off()
	
	cor(log10(rawExps+1))
	
	################################################################################
	#############  DE
	################################################################################
	
	
	DEGSeq.callDE(  exps+1 , condition=rep(c('DEL','control'),each=2) , countThres=0 ) 
	deRes = read.table("DEGseq_result/output_score.txt",sep="\t",header=T)
	deRes = deRes[ order(as.integer(as.character(deRes$GeneNames))) , ]
	info = data.frame(rawExps,exps,deRes[,c(4,7:8)])
	colnames(info) = c('flncDEL1','flncDEL2','flncCont1','flncCont2',
	'normalizedflncDEL1','normalizedflncDEL2','normalizedflncCont1','normalizedflncCont2',
	'log2FoldChange','pvalue','fdr')
	write.table(info,"DEGseq_result/output_info_shijia.txt",sep="\t",col.names=T,row.names=F,quote=F)
	
	
	
	################################################################################
	########## validate by adult
	################################################################################
	adultSamples= c("NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13")
	agffs = lapply( adultSamples , function(sample)
	{
	  files = dir( paste0(folder,"/",sample,"/",parai) )
	  gff = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
	  readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
	} )
	adultGff = unionGff( unionGff(agffs[[1]],agffs[[2]]) , agffs[[3]])
	
	validate = rep(0,length(allGff))
	validate[ matchGff(adultGff,allGff) ] = 1
	Name = reorder( paste0('isoform',1:length(allGff)) , de )
	validate2 = data.frame( Name , validate ) 
	heatmapMannualValidateColor(validate2)

	
	ids = sapply( translatedGffs , function(x) {
	  id = rep('NA',length(allGff))
  	id[ matchGff(x,allGff) ] = names(x)
    id
	})
	colnames(ids) = c("idDEL1","idDEL2","idCont1","idCont2")
	yy = data.frame(ids,info,adultValidate=validate)
	write.table(yy,"translatedGffs_id_exp_DE_adultvalidate.txt",col.names=T,row.names=F,sep="\t",quote=F)
}
