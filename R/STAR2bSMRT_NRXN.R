#' STAR2bSMRT_NRXN
#'
#' @param genomeDir 
#' @param genomeFasta 
#' @param phqv 
#' @param flnc 
#' @param nfl 
#' @param SR1 
#' @param SR2 
#' @param useSJout 
#' @param adjustNCjunc 
#' @param thresSR 
#' @param thresDis 
#' @param outputDir 
#' @param fixedMatchedLS 
#' @param fuzzyMatch
#' @param chrom 
#' @param s 
#' @param e 
#' @param cores 
#'
#' @return
#' @export
#'
#' @examples
#' 
STAR2bSMRT_NRXN <- function( genomeDir, genomeFasta, LRphqv=NULL, 
                        SR1, SR2=NULL, useSJout=TRUE,  adjustNCjunc=FALSE, 
                        thresSR, thresDis, outputDir, fixedMatchedLS=FALSE, fuzzyMatch=100, 
                        chrom=NULL , s=0 , e=Inf , cores=10 )
{

	library(Biostrings)
	library(foreach)
	library(doMC)
	registerDoMC(cores)
	
	############################################################################
	############   STARshort mapping and junction sites for short reads
	############################################################################
	
	SoutputDir = paste0(outputDir,"/SR")
	SRalignment = paste0(SoutputDir,"/alignments.bam")
	system( paste0( "mkdir -p " , SoutputDir ) )
	starShort( genomeDir , SR1 , SR2 , SoutputDir )
	
	if( useSJout )
	{
	  SRjunc = getJuncBySJout( SJout="SJ.out.tab", SoutputDir, chrom=chrom, s=s, e=e )
	} else {
	  SRjunc = getJuncBySam( SRalignment, SoutputDir, chrom=chrom, s=s, e=e )
	}
	
	
	############################################################################
	############   STARlong mapping and junction sites for long reads
	############################################################################
	
	############   for phqv	############   
	
	if( !is.null(LRphqv)  )
	{
	  LoutputDir = paste0(outputDir,"/LR")
	  LRalignment = paste0(LoutputDir,"/Aligned.out.sam")
	  system( paste0( "mkdir -p " , LoutputDir ) )
	  starLong( genomeDir=genomeDir , LR=LRphqv , outputDir=LoutputDir , cores=cores , SJ=NULL )
	  
	  LRread = getReadByJI( LRalignment , LoutputDir )
	  LRread = subset(LRread , start < 50150000 )
	  exp = phqvExp(LRphqv,LoutputDir)  # get coverage for all phqv 
	  LRread = merge( LRread , exp , by="id" ) 
	  LRread$coverage = LRread$full_length_coverage
	  LRinfo = getLRinfo( LRread ,  chrom=chrom , s=s , e=e )
	  LRread = LRinfo$LRread
	  LRjunc = LRinfo$LRjunc
	  LRtag = LRinfo$LRtag
	}

	
	############################################################################
	############   grid searching
	############################################################################
	
	if( fixedMatchedLS )
	{
	  matchedLS = matchLSjunc( LRjunc , SRjunc )
	} else {
	  matchedLS = NULL
	}
	
	score = gridSearch( LRjunc , SRjunc , thresSR , thresDis , adjustNCjunc , matchedLS , fuzzyMatch )
	
	ij = which( score==max(score) , arr.ind=T )
	ts = thresSR[ ij[1,1] ]
	td = thresDis[ ij[1,2] ]
	cat( ts , td , score[ij] , '\n ')
	
	correction = generateCorrectedIsoform( LRjunc , SRjunc, LRtag , LRread  , ts , td , matchedLS , fuzzyMatch )
	print(correction[[1]][c(2,3)])
	
	EoutputDir = paste0(outputDir,"/STAR2bSMRT")
	system( paste0( "mkdir -p " , EoutputDir ) )
	
	setwd( EoutputDir )
	genome = readDNAStringSet(genomeFasta)
	
	###############################################################################################################
	
	pdf( paste0( "JuncExp_LR_ts",ts,"_td",td,".pdf") )
	
	juncExp = do.call( rbind, lapply( correction , function(x) x$LSjuncCount ))
	lrCount = log10(juncExp$lrCount)
	srCount = log10(juncExp$srCount)
	juncCorr = cor.test(srCount,lrCount,method='spearman')$estimate
	cols = sapply( juncExp$motif , function(x) ifelse(x==0,1,2) )
	cols[juncExp$motif==1] = 3
	plot( lrCount , srCount , col=cols , pch=17 , main=paste0("JuncExp by Long and Short Reads: r=",signif(juncCorr,3)) ,  xlab="Log10 Long Read" , ylab="Log10 Short Read"  )
	abline(lm( srCount~lrCount ))
	
	par(mfrow=c(2,1))
	log10fc = lrCount - srCount
	JuncNames = paste(juncExp$start , juncExp$end)
	barplot( log10fc , cex.names=0.6 , col=cols , ylab="log10(lrCount/srCount)", names=JuncNames , las=3 )
	
	dev.off()
	
	###############################################################################################################
	gffName = paste0( "isoform_ts",ts,"_td",td,".gff")
	exonList = juncToExon( juncList=correction[[chrom]]$isoform , s=50149082 , e=51255411 , tag=correction[[chrom]]$normalizedIsoformCount )
	writeGff( isoform=exonList , file = gffName )
	
	###############################################################################################################
	seq = generateSeq( genome=genome , isoform=exonList )
	fastaName = paste0( "isoform_ts",ts,"_td",td,".fa")
	writeXStringSet( seq$dna , fastaName )
	#writeXStringSet( seq$dna[seq$translated] , fastaName )
	
	###############################################################################################################
	kallisto = kallistoQuant( fastaName , SR1 , SR2 , EoutputDir )
	
	Sexp = log10(kallisto$tpm+1)
	Lexp = log10(correction[['chr2']]$normalizedIsoformCount+1)
	LSQuantCorr = cor.test(Lexp,Sexp)$estimate
	LSQuantPval = cor.test(Lexp,Sexp)$p.val
	
	###############################################################################################################
	pdf( paste0( "Quant_LR_ts",ts,"_td",td,".pdf") )
	cols = sapply( seq$translated , function(x) ifelse(x,2,1) )
	plot( Lexp , Sexp , pch=16 , col=cols , main=paste0("Quantification by Long and Short Reads: r=",signif(LSQuantCorr,3)) ,  xlab="Log10 Long Read" , ylab="Log10 Short Read"  )
	abline(lm( Sexp~Lexp ))
	dev.off()
	
	###############################################################################################################
	pdf( "gridSeach.pdf" )
	heatmap( score , Rowv = NA, Colv = NA, scale='none' )
	dev.off()
	
	###############################################################################################################
	isoformNum = sum(sapply(correction,function(x)x$uniNum))
	isoformFrac = mean(sapply(correction,function(x)x$frac))
	info = data.frame( shortRead=ts , distance=td , isoformNum=isoformNum , isoformFrac=isoformFrac , translated=sum(seq$translated) , juncCorr , LSQuantCorr , LSQuantPval )
	write.table(info,"summary.txt",quote=F,sep="\t",col.names=T,row.names=F)
	

}



