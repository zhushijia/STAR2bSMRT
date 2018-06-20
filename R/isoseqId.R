#' isoseqId
#'
#' @param roi 
#' @param outputDir 
#'
#' @return
#' @export
#'
#' @examples
isoseqId = function( roi , outputDir )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/data/isoseqId.sh")
	output = paste0(outputDir,"/alignments.id")
	sh = paste( func , roi , output  )
	if( !file.exists(output) ) runSH(sh)

	read.table( output , sep="\t" , header=T )
}
