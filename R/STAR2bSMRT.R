#' STAR2bSMRT
#'
#' @param genomeDir 
#' @param LR 
#' @param SR1 
#' @param SR2 
#' @param outputDir 
#' @param adjustNCjunc 
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
STAR2bSMRT <- function( genomeDir , LR , SR1 , SR2 , outputDir , adjustNCjunc , chrom=NULL , s=0 , e=Inf , cores=1 )
{

	library(Biostrings)
	library(foreach)
	library(doMC)
	registerDoMC(cores)
	
	LoutputDir = paste0(outputDir,"/LR")
	SoutputDir = paste0(outputDir,"/SR")
	EoutputDir = paste0(outputDir,"/Exp")
	
	system( paste0( "mkdir -p " , LoutputDir ) )
	system( paste0( "mkdir -p " , SoutputDir ) )
	system( paste0( "mkdir -p " , EoutputDir ) )
	
	starShort( genomeDir , SR1 , SR2 , SoutputDir )
	starLong( genomeDir , LR , LoutputDir )
	genome = readDNAStringSet(ref)
	
	SRalignment = paste0(SoutputDir,"/alignments.bam")
	LRalignment = paste0(LoutputDir,"/Aligned.out.sam")
	
	SRjunc = getJunc( SRalignment , SoutputDir , chrom , s= 50147488 , e = 51259537 )
	LRinfo = getLRinfo( LRalignment , LR , LoutputDir , chrom , s= 50147488 , e = 51259537 )
	LRread = LRinfo$LRread
	LRjunc = LRinfo$LRjunc
	LRtag = LRinfo$LRtag
	
	#matchedLS = matchLSjunc( LRjunc , SRjunc )
	score = gridSearch( LRjunc , SRjunc , thresSR , thresDis , adjustNCjunc , matchedLS=NULL )
	
	ij = which( score==max(score) , arr.ind=T )
	ts = thresSR[ ij[1,1] ]
	td = thresDis[ ij[1,2] ]
	cat( ts , td , score[ij] , '\n ')
	
	correction = generateCorrectedIsoform( LRjunc , SRjunc, LRread  , ts , td , matchedLS=NULL )
	
	setwd( EoutputDir )
	
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
	writeGff( isoform=correction[[chrom]]$isoform , file = gffName , exp=correction[[chrom]]$exp , chrom='chr2' , s=50149082 , e=51255411 )
	
	###############################################################################################################
	seq = generateSeq( genome=genome , isoform=correction[[chrom]]$isoform , exp=correction[[chrom]]$exp , chrom='chr2' , s=50149082 , e=51255411  )
	fastaName = paste0( "isoform_ts",ts,"_td",td,".fa")
	writeXStringSet( seq$dna , fastaName )
	#writeXStringSet( seq$dna[seq$translated] , fastaName )
	
	###############################################################################################################
	kallisto = kallistoQuant( fastaName , SR1 , SR2 , EoutputDir )
	
	Sexp = log10(kallisto$tpm+1)
	Lexp = log10(correction[['chr2']]$exp+1)
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
	heatmap( P , Rowv = NA, Colv = NA, scale='none' )
	dev.off()

	###############################################################################################################
	isoformNum = sum(sapply(correction,function(x)x$num))
	isoformFrac = mean(sapply(correction,function(x)x$frac))
	info = data.frame( shortRead=ts , distance=td , isoformNum=isoformNum , isoformFrac=isoformFrac , translated=sum(seq$translated) , juncCorr , LSQuantCorr , LSQuantPval )
	write.table(info,"summary.txt",quote=F,sep="\t",col.names=T,row.names=F)
	

}



