#' kallistoQuant
#'
#' @param transcript 
#' @param SR1 
#' @param SR2 
#' @param outputDir 
#'
#' @return
#' @export
#'
#' @examples
kallistoQuant = function( transcript , SR1 , SR2 , outputDir )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/data/kallistoQuant.sh")
	sh = paste( func , transcript , SR1 , SR2 , outputDir  )
	runSH(sh)
	read.table( "output/abundance.tsv" , header=T )
}
