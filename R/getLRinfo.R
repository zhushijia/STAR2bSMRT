#' getLRinfo
#'
#' @param alignments 
#' @param phqv 
#' @param outputDir 
#' @param chrom 
#' @param s 
#' @param e 
#'
#' @return
#' @export
#'
#' @examples
getLRinfo = function( alignments , phqv=NULL , outputDir , chrom , s , e )
{
  LRread = getRead( alignments , outputDir )
  
  if( !is.null(phqv) )
  {
    exp = phqvExp( phqv , outputDir )
	  LRread = merge( LRread , exp , by="id" )
	  LRread$coverage = LRread$full_length_coverage # + LRread$non_full_length_coverage
  } else {
    LRread$coverage = 1
  }
  
	LRread = subset( LRread , junc!="-1" )

	if( !is.null(chrom) )
	{
		LRread = subset( LRread , chr %in% chrom & start>=s & start<50150000 & end<=e )
	  #LRread = subset( LRread , chr %in% chrom & start>=s & end<=e )
	}
	
	LRread = split( LRread , as.character(LRread$chr) )
	
	LRtag = lapply( LRread , function(read) 
			lapply( strsplit( as.character(read$junc),",") , function(x) { 
				i = 1:(length(x)/2); 
				paste( x[2*i-1] , x[2*i] ,sep="," ) 
		} ) )

	LRjunc = foreach( i = 1:length(LRread) ) %dopar%
	{
		read = LRread[[i]]
		tag = LRtag[[i]]
		len = sapply( tag, length )
		juncCount = rep( read$coverage , len )
		juncs = do.call(c,tag)
		uniJuncCount = tapply(juncCount,juncs,sum)
		start = sapply( strsplit(names(uniJuncCount),","), function(x) as.integer(x[1]) )
		end = sapply( strsplit(names(uniJuncCount),","), function(x) as.integer(x[2]) )
		data.frame(count=as.integer(uniJuncCount),chr='chr2',start=start,end=end,stringsAsFactors=F)
	}

	names(LRjunc) = names(LRread)

	list( LRread=LRread , LRjunc=LRjunc , LRtag=LRtag )
}


