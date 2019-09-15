#' matchLSjunc
#' match the junction site between long reads and short reads
#' 
#' @param LRjunc a list of data frame, indicating the splicing junction sites 
#' obtained from long reads
#' @param SRjunc a list of data frame, indicating the splicing junction sites 
#' obtained from short reads
#'
#' @return a list of data frame for the junction sites of each chromosome. The 
#' data frame indicates the locations of junction sites and the 
#' correspoding read count from both long and short reads
#' 
#' @export
#'
#' @examples
matchLSjunc = function( LRjunc , SRjunc )
{
	CHR = intersect( names(SRjunc) , names(LRjunc) )

	matched = foreach( k=1:length(CHR) ) %dopar%
	{
		chr = CHR[k]
		cat(chr,"\n")
		lrc = LRjunc[[chr]]
		src = SRjunc[[chr]]

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
	names(matched) = CHR

	matched

}

