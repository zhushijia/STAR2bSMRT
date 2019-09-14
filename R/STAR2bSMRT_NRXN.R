#' STAR2bSMRT_NRXN
#'
#'
#' @param genomeDir character value indicating the directory of STAR genome 
#' index for both STARlong and STARshort read mapping
#' @param genomeFasta character value indicating the fasta file of genome 
#' reference
#' @param phqv character value indicating the Isoseq polished high QV trascripts
#'  in fasta/fastq, where 
#' read counts for each transcript consensus should be saved in transcript names 
#' @param flnc character value indicating the Isoseq full-length non-chimeric
#'  reads in fasta/fastq format
#' @param nfl character value indicating the Isoseq non-full-length reads in 
#' fasta/fastq format
#' @param SR1 character value indicating the short read file in fastq format: 
#' single-end or paired-end R1
#' @param SR2 character value indicating the short read file in fastq format: 
#' paired-end R2
#' @param useSJout boolean value indicating whether to use the STARshort 
#' generated SJ.out.tab for splicing junction. If FALSE, STAR2bSMRT infer 
#' the splicing junction from bam files. By default, FALSE.
#' @param adjustNCjunc boolean value indicating whether to minimize the 
#' non-canonical junction sites. 
#' @param thresSR a vector of integers indicating the searching range for the 
#' number of short reads which support the splicing junction sites.
#' @param thresDis a vector of integers indicating the searching range for the 
#' tolerance distance between short read-derived splicing junction and long 
#' read-derived junction. STAR2bSMRT will correct the long read-derived 
#' junction to the short read-derived junction, if more short reads than 
#' defined thresSR support that short read-derived junction, and the distance 
#' between long and short read junctions is shorter than the defined thresDis.
#' @param outputDir character value indicating the direcotry where results are 
#' saved.
#' @param fixedMatchedLS boolean value indicating how often the distance is 
#' calculate betwen long read and short read-derived junction sites. If TRUE, 
#' only calculated once at the very beginning, which may save running time; 
#' otherwise, calculate repeatly after every long read correction. 
#' By default, FALSE.
#' @param fuzzyMatch integer value indicating the distance for fuzzyMatch
#' @param chrom character value indicating the chromosome of interest. By default, 
#' STAR2bSMRT works on the whole genome. 
#' @param s integeter value indicating the start position of the transcript of 
#' interest. This is useful for target Isoseq sequencing. 
#' @param e integeter value indicating the end position of the transcript of 
#' interest. This is useful for target Isoseq sequencing. 
#' @param cores integer value indicating the number of cores for parallel computing
#'
#' @return NULL
#' @export
#'
#' 
STAR2bSMRT_NRXN <- function( genomeDir, genomeFasta, LRphqv=NULL, LRflnc=NULL, LRnfl=NULL,
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
	tag = paste0( "exp", correction[[chrom]]$normalizedIsoformCount )
	exonList = juncToExon( juncList=correction[[chrom]]$isoform , s=50149082 , e=51255411 , tag=tag )

	###############################################################################################################
	seq = generateSeq( genome=genome , isoform=exonList )
	fastaName = paste0( "isoform_ts",ts,"_td",td,".fa")
	translated = sapply( seq$translated , function(x) ifelse( x , "translated" , "untranslated" )  )
	names(seq$rna) = paste(names(seq$rna),translated,sep="_")
	writeXStringSet( seq$rna , fastaName )
	
	gffName = paste0( "isoform_ts",ts,"_td",td,".gff")
	names(exonList) = paste( names(exonList) , translated , sep="_" )
	writeGff( isoform=exonList , file = gffName )
	
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



