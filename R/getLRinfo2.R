#' getLRinfo
#'
#' @param alignments 
#' @param phqv 
#' @param outputDir 
#' @param chrom 
#' @param s 
#' @param e 
#' @param jI 
#'
#' @return
#' @export
#'
#' @examples
getLRinfo <- function( alignments , phqv=NULL , outputDir , 
                      chrom=NULL , s=0 , e=Inf , jI=TRUE )
{
  library(reshape)
  
  if( jI )
  {
    LRread = getReadByJI( alignments , outputDir )
  } else {
    LRread = getReadByCigar( alignments , outputDir )
  }
  
  if( !is.null(phqv) )
  {
    exp = phqvExp( phqv , outputDir )
	  LRread = merge( LRread , exp , by="id" )
	  LRread$coverage = LRread$full_length_coverage 
	  # + LRread$non_full_length_coverage
  } else {
    LRread$coverage = 1
  }
  
	LRread = subset( LRread , junc!="-1" )
	LRread$smrtcell = sapply( strsplit(as.character(LRread$id),"/"), function(x) x[1] )
	uniSmrtcell = unique(LRread$smrtcell)
	
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
		juncCount = rep( read$coverage , len )
		smrtcell = rep( read$smrtcell , len )
		juncs = do.call(c,tag)
		
		info = data.frame(smrtcell,juncs,juncCount)
		pseudoInfo = data.frame( smrtcell=uniSmrtcell, 
		                         juncs= as.character(info$juncs[1]), juncCount=0 )
		info = rbind(pseudoInfo,info)
		
		count = cast( info , juncs ~ smrtcell, sum , value="juncCount")
		uniJunc = strsplit( as.character(count$juncs) , "," )
		start = sapply( uniJunc , function(x) as.integer(x[1]) )
		end = sapply( uniJunc , function(x) as.integer(x[2]) )
		data.frame( chr=names(LRread)[i], start=start, end=end, count[,-1],
		            stringsAsFactors=F)
	}

	names(LRjunc) = names(LRread)

	list( LRread=LRread , LRjunc=LRjunc , LRtag=LRtag )
	
}


