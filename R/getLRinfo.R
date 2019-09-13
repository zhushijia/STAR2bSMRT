#' getLRinfo
#'
#' @param LRread a list of data frame, indicating the information of long read,
#' comprising read id, chr, strand, start, end and junction sites
#' @param chrom a character value, incidating the chromosome of interest
#' @param s an integer value indicating the start site of gene of interest. This
#' is useful for targeted sequencing
#' @param e an integer value indicating the end site of gene of interest. This
#' is useful for targeted sequencing
#' 
#' @return a list of corrected isoforms, comprising the following items:
#' \itemize{
#'  \item {LRread} a list of data frame, indicating the information of long read,
#' comprising id, chr, strand, start, end, junc
#'  \item {LRjunc} a list of data frame, indicating the splicing junction sites 
#' obtained from long reads
#'  \item {LRtag} a list of character values, indicating the splicing junctions 
#' for each long read in the format of "junction_start,junction_end"
#'  
#' @export
#'
#' @examples
#' 
#' 
getLRinfo <- function( LRread, chrom=NULL, s=0, e=Inf )
{
  library(reshape)
  
  if( !"coverage" %in% colnames(LRread) )
  {
    LRread$coverage = 1
  }
    
  if( !"group" %in% colnames(LRread) )
  {
    LRread$group = "group1"
    uniGroup = "group1"
  } else {
    uniGroup = unique(LRread$group)
  }
  
	LRread = subset( LRread , junc!="-1" )
	
	if( !is.null(chrom) )
	{
		#LRread = subset( LRread , chr %in% chrom & start>=s & start<50150000 & end<=e )
	  LRread = subset( LRread , chr %in% chrom & start>=s & end<=e )
	}
	
	LRread = split( LRread , as.character(LRread$chr) )
	
	LRtag = lapply( LRread , function(read) 
			lapply( strsplit( as.character(read$junc),",") , function(x) { 
				i = 1:(length(x)/2); 
				paste( x[2*i-1] , x[2*i] ,sep="," ) 
		} ) )

	gc()
	
	LRjunc = foreach( i = 1:length(LRread) ) %dopar%
	{
	  cat(i,"\n")
		read = LRread[[i]]
		tag = LRtag[[i]]
		len = sapply( tag, length )
		juncCounts = rep( read$coverage , len )
		groups = rep( read$group , len )
		juncs = do.call(c,tag)
		
		info = data.frame(groups,juncs,juncCounts)
		pseudoInfo = data.frame( groups=uniGroup, 
		                         juncs= as.character(info$juncs[1]), juncCounts=0 )
		info = rbind(pseudoInfo,info)
		
		count = cast( info , juncs ~ groups, sum , value="juncCounts")
		uniJunc = strsplit( as.character(count$juncs) , "," )
		start = sapply( uniJunc , function(x) as.integer(x[1]) )
		end = sapply( uniJunc , function(x) as.integer(x[2]) )
		data.frame( chr=names(LRread)[i], start=start, end=end, count[,-1],
		            stringsAsFactors=F)
	}

	names(LRjunc) = names(LRread)

	list( LRread=LRread , LRjunc=LRjunc , LRtag=LRtag )
	
}


