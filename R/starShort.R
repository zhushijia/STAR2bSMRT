#' starShort
#'
#' @param genomeDir 
#' @param SR1 
#' @param SR2 
#' @param outputDir 
#' @param SJ 
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
