#' parseCigar
#'
#' @param cigar 
#'
#' @return
#' @export
#'
#' @examples
#' 
parseCigar = function( cigar )
{
	#cigar =  "7742M1I1954M1D1139M1I69M1I552M1I76M1I1964M1I4228M1D2023M1I2288M2D1326M1I1486M1I393M1I427M1I162M2D542M"
	cigar = as.character(cigar)
	pattern = "(\\d+M)|(\\d+I)|(\\d+D)|(\\d+N)|(\\d+S)|(\\d+H)|(\\d+P)|(\\d+=)|(\\d+X)"
	ind = gregexpr( pattern , cigar )[[1]]
	numL = as.numeric(ind)
	numR = numL + attr(ind, "match.length")- 2
	OpI  = numL + attr(ind, "match.length")- 1
	num  = mapply( function(s,e) as.integer(substr( cigar , s , e )) , numL , numR )
	Op   = sapply( OpI , function(s) substr( cigar , s , s ) )
	
	data.frame(num,Op)
	
}
