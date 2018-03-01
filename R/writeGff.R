#' writeGff
#'
#' @param isoform 
#' @param file 
#' @param exp 
#' @param chrom 
#' @param s 
#' @param e 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
writeGff <- function( isoform , file = "" )
{
  gff = list()
  for (i in 1:length(isoform))
  {
    exoni = isoform[[i]]
    len = nrow(exoni)
    chr = rep( as.character(exoni$chr[1]) , len+1 )
    generator = rep("SS",len+1)
    type = c( 'transcript' , rep( 'exon',len ) )
    start = c( exoni$start[1] , exoni$start )
    end = c( exoni$end[len] , exoni$end )
    strand = "-"
    id = paste0( 'gene_id "SS.1"; transcript_id "',names(isoform)[i],'";' )
    gff[[i]] = data.frame( chr,generator,type,start,end,".",strand,'.',id )
  }
  
  gff = do.call(rbind,gff)
  
  write.table( gff , file , col.names=F , row.names=F , sep="\t" , quote=F )
  
}

