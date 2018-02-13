#' matchGff
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
matchGff <- function( gff1 , gff2 )
{
  g1 = lapply( gff1,function(x) paste(paste( x[,1] , x[,2]  ),collapse="; ") )
  g2 = lapply( gff2,function(x) paste(paste( x[,1] , x[,2]  ),collapse="; ") )
  match(g1,g2)
}
