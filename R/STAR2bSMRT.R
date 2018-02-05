
#setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/R3/R")
#sapply( dir()[grep(".R",dir())] , source )



STAR2bSMRT = function( genomeDir , LR , SR1 , SR2 , outputDir , chrom=NULL , s=0 , e=Inf )
{
  
  
  STAR2bSMRT.dir ="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/R3"
  chrom = "chr2"
  s = 50147488
  e = 51259537
  cores = 30
  thresSR=c(10:30) 
  thresDis=c(1:20)
  
  
  
  
	library(foreach)
	library(doMC)
	registerDoMC(cores)
	
	LoutputDir = paste0(outputDir,"/LR")
	SoutputDir = paste0(outputDir,"/SR")

	system( paste0( "mkdir -p " , LoutputDir ) )
	system( paste0( "mkdir -p " , SoutputDir ) )

	# starLong( genomeDir , LR , LoutputDir )
	# starShort( genomeDir , SR1 , SR2 , SoutputDir )
	
	SRalignment = paste0(SoutputDir,"/alignments.bam")
	LRalignment = paste0(LoutputDir,"/Aligned.out.sam")
	
	SRjunc = getJunc( SRalignment , SoutputDir , chrom , s , e )
	LRinfo = getLRinfo( LRalignment , LR , LoutputDir , chrom , s , e )
	LRread = LRinfo$LRread
	LRjunc = LRinfo$LRjunc
	LRtag = LRinfo$LRtag
	
	matchedLS = matchLSjunc( LRjunc , SRjunc )
	
	P = gridSearch( LRjunc , SRjunc , matchedLS , thresSR=c(1:100) , thresDis=c(1:20) )
	ij = which( P==max(P) , arr.ind=T )
	ij
	
}






#LR.exp()
#gridSearch()
#LR.SR.JuncNum.ScatterPlot()
#LR.exp()
#generateFa()
#generateGff()
#mergeGff()

