#' juncToExon
#'
#' @param juncList 
#' @param s 
#' @param e 
#' @param tag
#'
#' @return
#' @export
#'
#' @examples
juncToExon = function( juncList , s , e , tag=c() )
{
	
  exons = lapply( juncList , function(junc) {
    chr = as.character(junc$chr[1])
    junc$start = junc$start-1
    junc$end = junc$end+1
    junc = rbind( data.frame(chr=chr,start=0,end=s) , junc, data.frame(chr=chr,start=e,end=Inf) )
    exon = data.frame( chr=junc$chr[-1] , start=junc$end[-nrow(junc)] , end=junc$start[-1] )
  } )
  
  names(exons) = paste0("S2bS",1:length(exons))
  
  if(length(tag)>0)
  {
    names(exons) = paste(names(exons),tag,sep="_")
  } 
  
  exons
}

