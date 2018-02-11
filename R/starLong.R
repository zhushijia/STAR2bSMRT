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
starLong = function( genomeDir , LR , outputDir )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/data/starLong.sh")
	sh = paste( func , genomeDir , LR , outputDir  )
	runSH(sh)
}
