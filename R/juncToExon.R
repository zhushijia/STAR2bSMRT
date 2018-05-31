#' juncToExon
#'
#' @param juncList 
#' @param s 
#' @param e 
#' @param exp 
#'
#' @return
#' @export
#'
#' @examples
juncToExon = function( juncList , s , e , exp=c() )
{
	
  exons = lapply( juncList , function(junc) {
    chr = as.character(junc$chr[1])
    junc$start = junc$start-1
    junc$end = junc$end+1
    junc = rbind( data.frame(chr=chr,start=0,end=s) , junc, data.frame(chr=chr,start=e,end=Inf) )
    exon = data.frame( chr=junc$chr[-1] , start=junc$end[-nrow(junc)] , end=junc$start[-1] )
  } )
  names(exons) = paste0("SS",1:length(exons),"_exp",exp)
  exons
}

