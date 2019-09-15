#' getReadByJI
#' get the long read information from JI tag in sam file
#' @param alignments character value indicating the location of alignment files,
#' such as .bam or .sam files
#' @param outputDir character value indicating the output directory for saving 
#' results
#'
#' @return a list of data frame, indicating the information of long read,
#' comprising id, chr, strand, start, end, junc
#' @export
#'
#' @examples
getReadByJI = function( alignments , outputDir )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/data/getRead.sh")
	output = paste0(outputDir,"/alignments.read")
	sh = paste( func , alignments , output  )
	if( !file.exists(output) ) runSH(sh)

	read.table( output , sep="\t" , header=T )
}

