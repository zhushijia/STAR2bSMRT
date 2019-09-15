#' phqvExp
#'
#' @param phqv character value indicating the Isoseq polished high QV trascripts
#'  in fasta/fastq, where 
#' @param outputDir character value indicating the direcotry where results are 
#' saved.
#'
#' @return data frame representing the read count for each phqv consensus
#' @export
#'
#' @examples
phqvExp = function( phqv , outputDir )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/data/phqvExp.sh")
	output = paste0(outputDir,"/phqv.exp")
	sh = paste( func , phqv , output  )
	if( !file.exists(output) ) runSH(sh)

	read.table( output , sep="\t" , header=T )
}
