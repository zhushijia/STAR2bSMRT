#' matchLSjuncOneChr
#'
#' @param lrc 
#' @param src 
#' @param fuzzyMatch 
#'
#' @return
#' @export
#'
#' @examples
#' 
matchLSjuncOneChr = function( lrc , src , fuzzyMatch=TRUE )
{
  
	SRmatch = matrix( ncol=2 , nrow=nrow(lrc) , data=0 )
	colnames(SRmatch) = c( "SRindex" , "LSdis" )

	for(i in 1:nrow(lrc))
	{
	  if( !fuzzyMatch )
	  {
	    dis = abs(src$start-lrc$start[i]) + abs(src$end-lrc$end[i])
	  } else {
	    Sdis = src$start-lrc$start[i]
	    Edis = src$end-lrc$end[i]
	    dis = abs(Edis-Sdis)
	  }
		
		ind = which.min(dis) 
		min = dis[ind]
		SRmatch[i,] = c(ind,min)
	}
	
	SRmatch
	
}

