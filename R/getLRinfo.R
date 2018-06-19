#' getLRinfo
#'
#' @param alignments a vector of character values representing the list of 
#' long read alignments
#' @param outputDir 
#' @param usejI a vector of boolean values with the same length of alignments 
#' representing whether read 
#' @param exp 
#' @param group a vector of character values representing either bin or batch. 
#' For Isoseq, it can represent the smrtcell ids, so as to perform normalization
#' across different smrtcells.
#' @param chrom 
#' @param s 
#' @param e 
#' 
#' @return
#' @export
#'
#' @examples
#' 
#' exp = phqvExp( phqv , outputDir )
#' exp = exp$full_length_coverage
#' group = sapply( strsplit(as.character(LRread$id),"/"), function(x) x[1] )
#' 
getLRinfo <- function( alignments, outputDir=dirname(alignments), 
                       usejI=rep(TRUE,length(alignments)), 
                       exp=as.list(rep(1,length(alignments))), 
                       group=as.list(rep("group1",length(alignments))), 
                       chrom=NULL, s=0, e=Inf )
{
  library(reshape)
  
  LRread <- lapply( seq_along(alignments) , function(i) {
    
    if( usejI[i] )
    {
      x = getReadByJI( alignments[i] , outputDir[i] )
    } else {
      x = getReadByCigar( alignments[i] , outputDir[i] )
    }
    
    x$coverage = exp[[i]]
    x$group = group[[i]]
    
    x
    
  })
  
  LRread = do.call(rbind, LRread)
	LRread = subset( LRread , junc!="-1" )
	uniGroup = unique(do.call(c,group))
	
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


