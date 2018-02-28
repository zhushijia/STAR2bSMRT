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
writeGff <- function( isoform , file = "" , exp , chrom , s , e )
{
  
  exons = lapply( isoform , function(junc) {
    chr = as.character(junc$chr[1])
    junc$start = junc$start-1
    junc$end = junc$end+1
    junc = rbind( data.frame(chr=chr,start=0,end=s) , junc, data.frame(chr=chr,start=e,end=Inf) )
    exon = data.frame( chr=junc$chr[-1] , start=junc$end[-nrow(junc)] , end=junc$start[-1] )
    exon
  } )
  names(exons) = paste0("SS",1:length(exons),"_exp",exp)
  
  gff = list()
  for (i in 1:length(exons))
  {
    exoni = exons[[i]]
    len = nrow(exoni)
    chr = rep( as.character(exoni$chr[1]) , len+1 )
    generator = rep("SS",len+1)
    type = c( 'transcript' , rep( 'exon',len ) )
    start = c( exoni$start[1] , exoni$start )
    end = c( exoni$end[len] , exoni$end )
    strand = "-"
    id = paste0( 'gene_id "SS.1"; transcript_id "',names(exons)[i],'";' )
    gff[[i]] = data.frame( chr,generator,type,start,end,".",strand,'.',id )
  }
  
  gff = do.call(rbind,gff)
  
  write.table( gff , file , col.names=F , row.names=F , sep="\t" , quote=F )
  
}

