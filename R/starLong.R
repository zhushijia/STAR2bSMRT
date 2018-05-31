#' starLong
#'
#' @param genomeDir 
#' @param LR 
#' @param outputDir 
#' @param cores 
#' @param SJ 
#'
#' @return
#' @export
#'
#' @examples
starLong = function( genomeDir , LR , outputDir , cores=10 , SJ=NULL )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	if( is.null(SJ) )
	{
	  func = paste0("source " , STAR2bSMRT.dir , "/data/starLong1.sh")
	  sh = paste( func , genomeDir , LR , outputDir , cores )
	  runSH(sh)
	} else {
	  func = paste0("source " , STAR2bSMRT.dir , "/data/starLong2.sh")
	  sh = paste( func , genomeDir , LR , outputDir , cores , SJ  )
	  runSH(sh)
	}
}
