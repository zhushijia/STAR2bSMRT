#' getJuncBySJout
#'
#' getJuncBySjout gets the unique junction sites and their numbers from the alignment files
#'  
#' @param SJout a data frame, indicating the SJ.out.tab file ouput by STAR mapping
#' @param outputDir a character value, indicating the directory for saving results
#' @param chrom a character value, incidating the chromosome of interest
#' @param s an integer value indicating the start site of gene of interest. This
#' is useful for targeted sequencing
#' @param e an integer value indicating the end site of gene of interest. This
#' is useful for targeted sequencing
#'
#' @return a list of data frame representing the splicing junction sites for 
#' each chromosome 
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
