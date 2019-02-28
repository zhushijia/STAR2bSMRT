library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")

plotPCA = function(gffs, exps)
{
  allGffs = NULL
  for(i in 1:length(gffs)) allGffs = unionGff(allGffs,gffs[[i]])
  
  binary = sapply( 1:length(gffs) , function(i) {
    flag = rep(0,length(allGffs))
    index = matchGff( allGffs, gffs[[i]] )
    flag[ !is.na(index) ] = 1
    flag
  })
  colnames(binary) = names(gffs)
  
  quantitative = sapply( 1:length(gffs) , function(i) {
    flag = rep(0,length(allGffs))
    index = matchGff( allGffs, gffs[[i]] )
    flag[ !is.na(index) ] = exps[[i]][ index[!is.na(index)] ]
    flag
  })
  colnames(quantitative) = names(gffs)
  
  pc = prcomp( t(binary) )$x
  plot( pc[,1] , pc[,2] , col="white",main="Isoform Yes or No" )
  text( pc[,1] , pc[,2] , colnames(binary) )
  dd <- dist( t(binary) , method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  plot(hc, main="Isoform Yes or No" )
  
  pc = prcomp( t(quantitative) )$x
  plot( pc[,1] , pc[,2] , col="white" ,main="Isoform Long Read Count" )
  text( pc[,1] , pc[,2] , colnames(quantitative) )
  dd <- dist( t(quantitative) , method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  plot(hc, main="Isoform Long Read Count" )
  
  pairs( log2(quantitative+1), pch=16, col='red'  )
  
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


heatmapMannualColor = function(X, case_col="darkred", cont_col="darkgrey", shared_col="darkorange" )
{
  library(ggplot2)
  library(reshape)
  library(plyr)
  X.m <- melt(X)
  X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
  gg <- ggplot(X.s, aes(x=variable, y=Name))
  gg <- gg + geom_tile(aes(fill = factor(value)), colour = "white")
  #gg <- gg + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred")
  gg <- gg + scale_fill_manual(values=c(case_col,"white",shared_col,cont_col))
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


drawSchematics = function( gffs, 
                           case_ind, cont_ind, valid_ind,
                           case_annotation, cont_annotation, 
                           case_col="darkred", cont_col="darkgrey", shared_col="darkorange",
                           threshold=7 , rankDE=FALSE, fileName,  outputDir )
{
    
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
  
	# normalized to sum
	rawExps = exps
	colnames(rawExps) = names(translatedGffs)
	
	normalizedFactor = colSums(exps)/min(colSums(exps[,c(case_ind,cont_ind)]))
	exps = sapply( 1:ncol(exps) , function(i) exps[,i]/normalizedFactor[i] )
	colnames(exps) = names(translatedGffs)
	exps = data.frame(exps)
	
	de = t( apply( exps , 1 , function(x) {
	    case = mean(x[case_ind])
	    cont = mean(x[cont_ind])
	    log2(case+1)-log2(cont+1)
	} ) )
	
	de = as.numeric(de)
	
	if (rankDE )
	{
	    Name = reorder( paste0('isoform',1:nrow(exps)) , de )
	} else {
	    Name = reorder( paste0('isoform',1:nrow(exps)) , rowSums(log10(exps[,c(case_ind,cont_ind)]+1)) )
	}
	
	logExps = data.frame( Name, log10(exps[,c(case_ind,cont_ind)]+1) )
	fracs = data.frame( Name , fracs ) 
	
	################################################################################
	########## validate by adult
	################################################################################
	
	if( length(valid_ind)>0 )
	{
	    validGff = NULL
	    for (i in valid_ind) validGff = unionGff(validGff,translatedGffs[[i]])
	    validate = rep(0,length(allGff))
	    validate[ matchGff(validGff,allGff) ] = 1
	    validate2 = data.frame( Name , validate ) 
	    validateName = paste( c("validateBy", names(gffs)[valid_ind]), collapse="_")
	}
	
	################################################################################
	#############  isoform schematics
	################################################################################
	system( paste( "mkdir -p" , outputDir ) )
	setwd(outputDir)
	
	pdfName = paste( c(fileName,
	                   case_col,length(caseUniGff),"for",names(gffs)[case_ind],
	                   cont_col,length(contUniGff),"for",names(gffs)[cont_ind],
	                   shared_col,length(commonGff),"forShared",
                     ifelse( length(valid_ind)>0,validateName,""),
                     ifelse(rankDE,"rankByDEG","rankByExp"), 
                     "_threshold", threshold,".pdf"), collapse="_" )

	pdf( pdfName ) #,h=10,w=4)
	
	heatmapMannualColor(fracs, case_col, cont_col, shared_col)
	heatmapGradientColor(logExps)
	
	########## validate by others
	if( length(valid_ind)>0 ) heatmapMannualValidateColor(validate2)
	
	cased = matchGff(caseUniGff,allGff)
	contd = matchGff(contUniGff,allGff)
	shared = matchGff(commonGff,allGff)
	
	cols = rep('',nrow(exps))
	cols[cased] = case_col
	cols[contd] = cont_col
	cols[shared] = shared_col
	
	barplot(de[order(Name)],col=cols[order(Name)],horiz=T,border='white')
	
	plotPCA(translatedGffs, translatedExps)
	
	dev.off()
	
	#casecontGff = list(caseGff, contGff)
	#names(casecontGff) = c( paste( names(gffs)[case_ind], collapse=" " ),
	#                        paste( names(gffs)[cont_ind], collapse=" " ) )
	#vennDiagramGff( casecontGff, gsub('pdf','png',pdfName) )
	
}


################################################################################
#############  Fig 2. 2 hiPSC control; 3 Fetal; 3 Adults
################################################################################

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples= c("581","641","2607","553","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks") #
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Cont3","Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") #
case_ind = c(3,4)
cont_ind = c(6:11)
valid_ind = NULL
case_col="darkgrey"
cont_col="darkgreen"
shared_col="darkorange"
human_annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
case_annotation = human_annotation
cont_annotation = human_annotation
threshold=7
rankDE=FALSE
fileName = "IsoformSchematics_Fig2_2Cont_3Fetal_3Adult"
outputDir = paste0(folder,"/comparisonResult","/",parai)

gffs = lapply( samples , function(sample)
{
  files = dir( paste0(folder,"/",sample,"/",parai) )
  gff = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
  readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(gffs) = sampleNames

drawSchematics( gffs, case_ind, cont_ind, valid_ind, 
                case_annotation, cont_annotation, 
                case_col, cont_col, shared_col,
                threshold, rankDE, fileName, outputDir )


################################################################################
#############  Fig 3 2 case vs 2 cont
################################################################################

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
samples= c("581","641","2607","553","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13")
sampleNames = c("p3Del1","p3Del2","Cont1","Cont2","Adult1","Adult2","Adult3")
case_ind = c(1,2)
cont_ind = c(3,4)
valid_ind = c(5:7)
case_col="darkred"
cont_col="darkgrey"
shared_col="darkorange"
human_annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
case_annotation = human_annotation
cont_annotation = human_annotation
threshold=7
rankDE=TRUE
fileName = "IsoformSchematics_Fig3_2Case_2Cont"
outputDir = paste0(folder,"/comparisonResult","/",parai)

gffs = lapply( samples , function(sample)
{
  files = dir( paste0(folder,"/",sample,"/",parai) )
  gff = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
  readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(gffs) = sampleNames

drawSchematics( gffs, case_ind, cont_ind, valid_ind, 
                case_annotation, cont_annotation, 
                case_col, cont_col, shared_col,
                threshold, rankDE, fileName, outputDir )


################################################################################
#############  Fig 4. 553 GABA GLUT
################################################################################

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples = c("553_GABA","553_NGN2","553")
sampleNames = c("553_GABA","553_NGN2","553_FBN")
case_ind = c(1)
cont_ind = c(2)
valid_ind = c(3)
case_col="darkred"
cont_col="darkgrey"
shared_col="darkorange"
human_annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
case_annotation = human_annotation
cont_annotation = human_annotation
threshold=7
rankDE=FALSE
fileName = "IsoformSchematics_Fig4_553NGN2_GABA_FBN"
outputDir = paste0(folder,"/comparisonResult","/",parai)

gffs = lapply( samples , function(sample)
{
  files = dir( paste0(folder,"/",sample,"/",parai) )
  gff = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
  readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(gffs) = sampleNames

drawSchematics( gffs, case_ind, cont_ind, valid_ind, 
                case_annotation, cont_annotation, 
                case_col, cont_col, shared_col,
                threshold, rankDE, fileName, outputDir )

	


################################################################################
#############  Supple 3 Fetal; 3 Adults
################################################################################

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
samples= c("NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal","fetal_23wks","fetal_3wks") #
sampleNames = c("Adult1","Adult2","Adult3","Fetal1","Fetal2","Fetal3") #
case_ind = c(1:3)
cont_ind = c(4:6)
valid_ind = NULL
case_col="darkred"
cont_col="darkgrey"
shared_col="darkorange"
human_annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
case_annotation = human_annotation
cont_annotation = human_annotation
threshold=7
rankDE=FALSE
fileName = "IsoformSchematics_Supp_3Fetal_3Adult"
outputDir = paste0(folder,"/comparisonResult","/",parai)

gffs = lapply( samples , function(sample)
{
  files = dir( paste0(folder,"/",sample,"/",parai) )
  gff = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
  readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(gffs) = sampleNames

drawSchematics( gffs, case_ind, cont_ind, valid_ind, 
                case_annotation, cont_annotation, 
                case_col, cont_col, shared_col,
                threshold, rankDE, fileName, outputDir )



