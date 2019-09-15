#' salmonQuant
#' perform the salmon transcript quantification
#' @param transcript character value indicating the transcriptome file
#' @param SR1 character value indicating the short read file in fastq format: 
#' single-end or paired-end R1
#' @param SR2 character value indicating the short read file in fastq format: 
#' paired-end R2
#' @param outputDir character value indicating the direcotry where results are 
#' saved.
#'
#' @return data frame representing the expression for each transcript
#' @export
#'
#' @examples
salmonQuant = function( transcript , SR1 , SR2 , outputDir )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/data/kallistoQuant.sh")
	sh = paste( func , transcript , SR1 , SR2 , outputDir  )
	runSH(sh)
	read.table( "output/quant.sf" , header=T )
}
