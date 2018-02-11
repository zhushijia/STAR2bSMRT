#' generateCorrectedIsoform
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
generateCorrectedIsoform = function( LRjunc , SRjunc, LRread , matchedLS , ts , td )
{
	
	CHR = intersect( names(SRjunc) , names(LRjunc) )

	allCorrectIsoform = foreach( k=1:length(CHR) ) %dopar% {
		
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
		
		correctIsoformExp = tapply( read$full_length_coverage[range] , isoform , sum )
		srindex = lapply( strsplit( names(correctIsoformExp) , "_" ) , function(x) as.integer(x) )
		correctIsoform = lapply( srindex , function(ind) { 
			x=src[ ind , c('chr','start','end') ]
			x[ order(x$start), ]
			} )
		
		correctIsoformNum = length(unique(isoform))
		correctIsoformFrac = sum(correctIsoformExp)/sum(read$full_length_coverage) 
		
		list( isoform=correctIsoform , num=correctIsoformNum , exp=as.numeric(correctIsoformExp) , frac=correctIsoformFrac )
		
	}

	names(allCorrectIsoform) = CHR
  
	allCorrectIsoform
	
}
