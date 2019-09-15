#' matchLSjuncOneChr
#'
#' @param lrc 
#' @param src 
#' @param fuzzyMatch integer value indicating the distance for fuzzymatching 
#'
#' @return
#' @export
#'
#' @examples
#' 
matchLSjuncOneChr = function( lrc , src , fuzzyMatch=100 )
{
  
	SRmatch = matrix( ncol=2 , nrow=nrow(lrc) , data=0 )
	colnames(SRmatch) = c( "SRindex" , "LSdis" )

	for(i in 1:nrow(lrc))
	{
	  
    Sdis = src$start-lrc$start[i]
    Edis = src$end-lrc$end[i]
    dis1 = abs(src$start-lrc$start[i]) + abs(src$end-lrc$end[i])
    dis2 = abs(Edis-Sdis)
    
    mindis1 = min(  dis1  )
    inddis1 = which.min(dis1) 
    
    fuzzyRange = ( abs(Sdis)<fuzzyMatch & abs(Edis)<fuzzyMatch )
    if( sum(fuzzyRange)>0 )
    {
      mindis2=min( dis2[ fuzzyRange ])
      inddis2=which( fuzzyRange & dis2==mindis2 )[1]
    } else {
      mindis2=Inf
      inddis2=NA
    }
    
    if( mindis1<=mindis2 )
    {
      ind = inddis1
      min = mindis1
    } else {
      ind = inddis2
      min = mindis2
    }
    
		SRmatch[i,] = c(ind,min)
	}
	
	SRmatch
	
}

