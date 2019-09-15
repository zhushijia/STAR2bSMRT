#' unionGff
#' unions two isoform gff files
#' 
#' @param gff1 a list of data frame representing the junctions for each isoform
#' @param gff2 a list of data frame representing the junctions for each isoform
#'
#' @return a combined list of unique isoforms
#' @export
#'
#' @examples
#' 
#' 
unionGff <- function( gff1 , gff2 )
{
  index = matchGff(gff1 , gff2)
  index = index[!is.na(index)]
  index2 = setdiff( seq_along(gff2) , index  )
  c( gff1 , gff2[index2] ) 
}
