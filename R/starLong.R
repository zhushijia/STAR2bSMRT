#' starLong
#' for starlong mapping
#' @param genomeDir character value indicating the location of STAR genome indexes
#' @param LR character value indicating the long read file
#' @param outputDir character value indicating the directory of output
#' @param cores integer value indicating the number of cores for parallel computing
#' @param SJ character value for 2pass STAR mapping
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
