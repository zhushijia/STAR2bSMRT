#' matchLSjuncOneChr
#'
#' @param lrc 
#' @param src 
#'
#' @return
#' @export
#'
#' @examples
#' 
matchLSjuncOneChr = function( lrc , src )
{
  
	SRmatch = matrix( ncol=2 , nrow=nrow(lrc) , data=0 )
	colnames(SRmatch) = c( "SRindex" , "LSdis" )

	for(i in 1:nrow(lrc))
	{
		dis = abs(src$start-lrc$start[i]) + abs(src$end-lrc$end[i])
		ind = which.min(dis) 
		min = dis[ind]
		SRmatch[i,] = c(ind,min)
	}
	
	SRmatch
	
}

