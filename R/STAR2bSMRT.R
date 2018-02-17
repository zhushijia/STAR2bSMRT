
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
	
	matchedLS = matchLSjunc( LRjunc , SRjunc )
	P = gridSearch( LRjunc , SRjunc , matchedLS , thresSR , thresDis )
	
	ij = which( P==max(P) , arr.ind=T )
	ts = thresSR[ ij[1,1] ]
	td = thresDis[ ij[1,2] ]
	cat( ts , td , P[ij] , '\n ')
	
	correction = generateCorrectedIsoform( LRjunc , SRjunc, LRread , matchedLS , ts , td )
	
	setwd( EoutputDir )
	
	###############################################################################################################
	
	pdf( paste0( "JuncExp_LR_ts",ts,"_td",td,".pdf") )
	
	juncExp = do.call( rbind, lapply( correction , function(x) x$LRjuncCount ))
	lrCount = log10(juncExp$lrCount)
	srCount = log10(juncExp$srCount)
	juncCorr = cor.test(srCount,lrCount)$estimate
	cols = sapply( juncExp$motif , function(x) ifelse(x==0,1,2) )
	cols[juncExp$motif==1] = 3
	plot( lrCount , srCount , col=cols , pch=16 , main=paste0("JuncExp by Long and Short Reads: r=",signif(juncCorr,3)) ,  xlab="Log10 Long Read" , ylab="Log10 Short Read"  )
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
	seq = generateSeq( genome , isoform=correction[[chrom]]$isoform , exp=correction[[chrom]]$exp , chrom='chr2' , s=50149082 , e=51255411 )
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
	info = data.frame( shortRead=ts , distance=td , isoformNum=isoformNum , isoformFrac=isoformFrac , juncCorr , LSQuantCorr , LSQuantPval )
	write.table(info,"summary.txt",quote=F,sep="\t",col.names=T,row.names=F)
	
	

}




#setwd(SoutputDir)
#system( paste0( "samtools view alignments.bam " , chrom,":",s,"-",e," > cut.bam" ) )
#system( "samtools sort -n cut.bam cut.bam.qsort" )
#system( "bedtools bamtofastq -i cut.bam.qsort.bam -fq cut.R1.fastq -fq2 cut.R2.fastq")





