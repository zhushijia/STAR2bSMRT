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
  g1 = lapply( gff1,function(x) paste( x$chr[1] , paste( paste( x$start , x$end  ),collapse="; "),sep=": " ) )
  g2 = lapply( gff2,function(x) paste( x$chr[1] , paste( paste( x$start , x$end  ),collapse="; "),sep=": " ) )
  match(g1,g2)
}
