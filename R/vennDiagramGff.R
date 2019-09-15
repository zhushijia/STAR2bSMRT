#' vennDiagramGff
#' draw the venn diagram of a list of gffs
#' @param gffList a list of gffs representing gffs for different samples
#' @param filename character value representing the file name of VennDiagram
#'
#' @return
#' @export
#'
#' @examples
#' 
vennDiagramGff <- function(gffList,filename)
{
  library(VennDiagram)
  toTag = function(x) paste( x$chr[1] , paste( paste( x$start , x$end  ),collapse="; "),sep=": " )
  tagList = lapply( gffList , function(y) sapply( y , toTag ) )
  tagAll = unique( do.call(c, tagList ))
  indexList <- lapply( tagList , function(x) which(tagAll%in%x) )
  names(indexList) = names(gffList)
  venn.diagram(indexList, filename=filename , height = 2000, width = 2000, fill=1:length(gffList) , alpha=0.4 ,imagetype = "png")
}
