#' starShort
#' for STARshort mapping
#' @param genomeDir character value indicating the location of STAR genome indexes 
#' @param SR1 character value indicating the short read file in fastq format: 
#' single-end or paired-end R1
#' @param SR2 character value indicating the short read file in fastq format: 
#' paired-end R2
#' @param outputDir character value indicating the directory of output
#' @param SJ character value for 2pass STAR mapping
#'
#' @return
#' @export
#'
#' @examples
starShort = function( genomeDir , SR1 , SR2 , outputDir , SJ=NULL )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	if( is.null(SJ) )
	{
	  func = paste0("source " , STAR2bSMRT.dir , "/data/starShort1.sh")
	  sh = paste( func , genomeDir , SR1 , SR2 , outputDir  )
	  runSH(sh)
	} else {
	  func = paste0("source " , STAR2bSMRT.dir , "/data/starShort2.sh")
	  sh = paste( func , genomeDir , SR1 , SR2 , outputDir , SJ )
	  runSH(sh)
	}
}
