#' readGff
#'
#' @param file 
#' @param chrom 
#' @param s 
#' @param e 
#'
#' @return
#' @export
#'
#' @examples
readGff <- function( file , chrom , s , e )
{
  
  gff = read.table(file,sep='\t')
  colnames(gff)[c(1,3:5,7)] = c('chr','type','start','end','strand')
  exon = subset( gff , type=='exon' & start>=s & end<=e )
  
  if( !is.null(chrom) )
    exon = subset( exon , chr==chrom )
  
  split( exon[,c(1,4:5)] , exon[,9] )
  
}
