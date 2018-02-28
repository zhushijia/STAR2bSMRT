#' unionGff
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
  c( gff1 , gff2[-index] ) 
}
