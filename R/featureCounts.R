#' featureCounts
#'
#' @param gtf 
#' @param bam 
#' @param outputDir 
#' @param feature 
#' @param cores 
#'
#' @return
#' @export
#'
#' @examples
#' 
featureCounts <- function( gtf , bam , outputDir , feature='exon' , cores=10 )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	
	for( fi in feature )
	{
		func = paste0("source " , STAR2bSMRT.dir , "/data/featureCounts.sh")
		sh = paste( func , gtf , bam , outputDir , fi , cores )
		runSH(sh)
	}
	
}

