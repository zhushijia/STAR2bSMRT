
#setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/R3/R")
#sapply( dir()[grep(".R",dir())] , source )

#library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone/setup")

STAR2bSMRT = function( genomeDir , LR , SR1 , SR2 , outputDir , chrom=NULL , s=0 , e=Inf )
{
  
  
  ref = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
  chrom = "chr2"
  s = 50147488
  e = 51259537
  cores = 30
  thresSR=c(1:100) 
  thresDis=c(1:30)
  
  
  
  
  
	library(foreach)
	library(doMC)
	registerDoMC(cores)
	library(Biostrings)
	
	LoutputDir = paste0(outputDir,"/LR")
	SoutputDir = paste0(outputDir,"/SR")
	EoutputDir = paste0(outputDir,"/Exp")
	
	system( paste0( "mkdir -p " , LoutputDir ) )
	system( paste0( "mkdir -p " , SoutputDir ) )
	system( paste0( "mkdir -p " , EoutputDir ) )
	
	#starLong( genomeDir , LR , LoutputDir )
	#starShort( genomeDir , SR1 , SR2 , SoutputDir )
	
	SRalignment = paste0(SoutputDir,"/alignments.bam")
	LRalignment = paste0(LoutputDir,"/Aligned.out.sam")
	
	SRjunc = getJunc( SRalignment , SoutputDir , chrom , s= 50147488 , e = 51259537 )
	LRinfo = getLRinfo( LRalignment , LR , LoutputDir , chrom , s= 50147488 , e = 51259537 )
	LRread = LRinfo$LRread
	LRjunc = LRinfo$LRjunc
	LRtag = LRinfo$LRtag
	
	matchedLS = matchLSjunc( LRjunc , SRjunc )
	P = gridSearch( LRjunc , SRjunc , matchedLS , thresSR , thresDis )
	
	ij = which( P==max(P) , arr.ind=T )
	ts = thresSR[ ij[1,1] ]
	td = thresDis[ ij[1,2] ]
	cat( ts , td , P[ij] , '\n ')
	
	correction = generateCorrectedIsoform( LRjunc , SRjunc, LRread , matchedLS , ts , td )
	table(sapply(correction[[1]]$isoform,function(x)x[1,2]))
	table(sapply(correction[[1]]$isoform,function(x)x[nrow(x),3]))
	# genome = readDNAStringSet(ref)
	
	setwd( EoutputDir )
	seq = generateSeq( genome , isoform=correction[[chrom]]$isoform , exp=correction[[chrom]]$exp , chrom='chr2' , s=50149082 , e=51255411 )
	fastaName = paste0( "isoform_ts",ts,"_td",td,"fa")
	writeXStringSet( seq$dna , fastaName )
	#writeXStringSet( seq$dna[seq$translated] , "isoform.fa" )
	
	kallisto = kallistoQuant( fastaName , SR1 , SR2 , EoutputDir )
	
	Sexp = log10(kallisto$tpm+1)
	Lexp = log10(correction[['chr2']]$exp+1)
	cor.test(Sexp,Lexp)
	
	pdf( paste0( "Quant_LR_ts",ts,"_td",td,".pdf") )
	cols = sapply( seq$translated , function(x) ifelse(x,2,1) )
	plot( Lexp , Sexp , pch=16 , col=cols , main="Quantification by Long and Short Reads" ,  xlab="Log10 Long Read" , ylab="Log10 Short Read"  )
	abline(lm( Sexp~Lexp ))
	dev.off()
	
	pdf( "gridSeach.pdf" )
	heatmap( P , Rowv = NA, Colv = NA, scale='none' )
	dev.off()
	

}






#LR.exp()
#gridSearch()
#LR.SR.JuncNum.ScatterPlot()
#LR.exp()
#generateFa()
#generateGff()
#mergeGff()

