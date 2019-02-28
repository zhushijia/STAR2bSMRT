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
  X.m <- melt(X)
  X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
  gg <- ggplot(X.s, aes(x=variable, y=Name))
  gg <- gg + geom_tile(aes(fill = factor(value)), colour = "white")
  #gg <- gg + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred")
  gg <- gg + scale_fill_manual(values=c("darkred","white","darkorange","darkgrey"))
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




annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)

library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"

drawSchematics = function( samples, sampleNames, case_ind, cont_ind, valid_ind, fileName , threshold , rankDE=FALSE )
{
    
    parai="adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
    
	cat(parai,"\n")
    
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
	
	caseGff = NULL
	for (i in case_ind) caseGff = unionGff(caseGff,translatedGffs[[i]])
	contGff = NULL
	for (i in cont_ind) contGff = unionGff(contGff,translatedGffs[[i]])
	
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
	
	de = t( apply( exps , 1 , function(x) {
	    case = mean(x[case_ind])
	    cont = mean(x[cont_ind])
	    fc = log2(case+1)-log2(cont+1)
	} ) )
	
	de = as.numeric(de)
	
	
	if (rankDE )
	{
	    Name = reorder( paste0('isoform',1:nrow(exps)) , de )
	} else {
	    Name = reorder( paste0('isoform',1:nrow(exps)) , rowSums(log10(exps[,c(case_ind,cont_ind)]+1)) )
	}
	
	logExps = data.frame( Name, log10(exps+1) )
	fracs = data.frame( Name , fracs ) 
	
	outputDir = paste0(folder,"/comparisonResult","/",parai)
	system( paste( "mkdir -p" , outputDir ) )
	setwd(outputDir)
    
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
	}
	
	################################################################################
	#############  isoform schematics
	################################################################################
	
	pdf(  paste0(fileName,"_threshold",threshold,".pdf") ) #,h=10,w=4)
	
	
	heatmapMannualColor(fracs)
	heatmapGradientColor(logExps)
	########## validate by others
	if( length(valid_ind)>0 ) heatmapMannualValidateColor(validate2)
	
	cased = matchGff(caseUniGff,allGff)
	contd = matchGff(contUniGff,allGff)
	shared = matchGff(commonGff,allGff)
	cols = rep('',nrow(exps))
	cols[shared] = "darkorange"
	cols[cased] = "darkred"
	cols[contd] = "darkgrey"
	
	barplot(de[order(Name)],col=cols[order(Name)],horiz=T,border='white')
	
	dev.off()
	
	
	

}
	
drawSchematics( samples = c("553_GABA","553_NGN2","553") ,
                sampleNames = c("553_GABA","553_NGN2","553_FBN"), 
                case_ind = c(1) , cont_ind = c(2) , valid_ind = c(3) , 
                fileName="GABA_GLUT_FBN_553_schematics", threshold=7, rankDE=FALSE )

	
