#' matchGff
#' finds the correspondence between two gffs
#' 
#' @param gff1 a list of data frame representing the junctions for each isoform
#' @param gff2 a list of data frame representing the junctions for each isoform
#'
#' @return a vector of integer values indicating the index of gff2 matching gff1
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
