#' stringTie
#'
#' @param alignment 
#' @param gff 
#' @param outputDir 
#'
#' @return
#' @export
#'
#' @examples
#' 
stringTie = function( alignments , gff , outputDir )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/data/stringTie.sh")
	sh = paste( func , alignments , gff , outputDir  )
	runSH(sh)
}
