#' starLong
#'
#' @param genomeDir 
#' @param LR 
#' @param outputDir 
#'
#' @return
#' @export
#'
#' @examples
#' 
starLong = function( genomeDir , LR , outputDir , SJ=NULL )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	if( is.null(SJ) )
	{
	  func = paste0("source " , STAR2bSMRT.dir , "/data/starLong1.sh")
	  sh = paste( func , genomeDir , LR , outputDir  )
	  runSH(sh)
	} else {
	  func = paste0("source " , STAR2bSMRT.dir , "/data/starLong2.sh")
	  sh = paste( func , genomeDir , LR , outputDir , SJ )
	  runSH(sh)
	}
}
