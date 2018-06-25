#' getJuncBySJout
#'
#' getJuncBySjout gets the unique junction sites and their numbers from the alignment files
#'  
#' @param SJout a character value, indicating the SJ.out.tab file 
#' @param outputDir a character value, indicating the directory whether to save the results
#' @param chrom a vector of character values, incidating the chromosomes of insterest
#' @param s an integer value indicating the start site
#' @param e an integer value indicating the end site
#'
#' @return
#' @export
#'
#' @examples
#' 
getJuncBySJout = function( SJout="SJ.out.tab" , outputDir , chrom=NULL ,  s=0 , e=Inf )
{
  SJ.out.tab = read.table( paste0(outputDir,"/",SJout) , sep="\t")
  junc = with( SJ.out.tab , data.frame(count=V7, chr=V1, start=V2, end=V3, motif=V5) )
  
  if( !is.null(chrom) )
  {
    junc = subset( junc , chr %in% chrom & start>=s & end<=e )
  }
  
  split( junc , as.character(junc$chr) )
  
}
