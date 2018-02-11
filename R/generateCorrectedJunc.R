#' generateCorrectedJunc
#'
#' @param LRjunc 
#' @param SRjunc 
#' @param LRread 
#' @param matchedLS 
#' @param ts 
#' @param td 
#'
#' @return
#' @export
#'
#' @examples
generateCorrectedJunc = function( LRjunc , SRjunc, LRread , matchedLS , ts , td )
{
	
	CHR = intersect( names(SRjunc) , names(LRjunc) )

	Junc = foreach( k=1:length(CHR) ) %dopar% {
		
		chr = CHR[k]
		lrc = LRjunc[[chr]]
		src = SRjunc[[chr]]
		tag = LRtag[[chr]]
		read = LRread[[chr]]
		SRmatch = matchedLS[[chr]]

		index = which( SRmatch[,'LSdis']<=td & src$count[SRmatch[,'SRindex']]>=ts ) 
		LRcorres = data.frame(lrc[index,],SRmatch[index,])
		correctTag = paste(LRcorres$start , LRcorres$end , sep="," )
		range = sapply( tag , function(x) all(x%in%correctTag) )
		isoform = sapply( tag[range] , function(x) paste( sort(subset(LRcorres,correctTag%in%x)$SRindex) , collapse="_" ) )
		correctExp = tapply( read$full_length_coverage[range] , isoform , sum )
		srindex = lapply( strsplit( names(correctExp) , "_" ) , function(x) as.integer(x) )
		correctJunc = lapply( srindex , function(ind) { 
			x=src[ ind , c('chr','start','end') ]
			x[ order(x$start), ]
			} )
		correctJunc
	}

	names(Junc) = CHR

	length(unique(isoform))
	sum(range)
	sum(read$full_length_coverage[range])/sum(read$full_length_coverage)

}
