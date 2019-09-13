#' getJuncBySam
#'
#' getJuncBySam gets the unique junction sites and their numbers from the alignment files
#'  
#' @param alignments a character value, indicating the alignment file, .sam for .bam 
#' @param outputDir a character value, indicating the directory whether to save the results
#' @param chrom a vector of character values, incidating the chromosomes of insterest
#' @param s an integer value indicating the start site
#' @param e an integer value indicating the end site
#'
#' @return a list of data frame representing the splicing junction sites for 
#' each chromosome 
#' @export
#'
#' @examples
getJuncBySam = function( alignments , outputDir , chrom=NULL ,  s=0 , e=Inf )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/data/getJunc.sh")
	output = paste0( outputDir , "/alignments.junc" )
	sh = paste( func , alignments , output  )
	if( !file.exists(output) ) runSH(sh)

	junc = read.table( output , sep="\t" , header=T )

	if( !is.null(chrom) )
	{
		junc = subset( junc , chr %in% chrom & start>=s & end<=e )
	}
	
	split( junc , as.character(junc$chr) )

}
