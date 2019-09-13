#' unionGff
#' unions two isoform gff files
#' 
#' @param gff1 
#' @param gff2 
#'
#' @return
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
