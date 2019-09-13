#' readGff
#' for reading isoform gff file
#' 
#' @param file character value indicating the name of the gff file
#' @param chrom character value indicating the chromosome of interest. By default, 
#' STAR2bSMRT works on the whole genome. 
#' @param s integeter value indicating the start position of the transcript of 
#' interest. This is useful for target Isoseq sequencing. 
#' @param e integeter value indicating the end position of the transcript of 
#' interest. This is useful for target Isoseq sequencing. 
#'
#' @return a list of isoforms 
#' @export
#'
#' @examples
readGff <- function( file , chrom=NULL , s=0 , e=Inf )
{
  
  gff = read.table(file,sep='\t')
  colnames(gff)[c(1,3:5,7)] = c('chr','type','start','end','strand')
  exon = subset( gff , type=='exon' & start>=s & end<=e )
  
  if( !is.null(chrom) )
    exon = subset( exon , chr==chrom )
  
  split( exon[,c(1,4:5)] , exon[,9] )
  
}
