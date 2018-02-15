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
		#LRcorres = subset(LRcorres,count>3)
		correctTag = paste(LRcorres$start , LRcorres$end , sep="," )
		range = sapply( tag , function(x) all(x%in%correctTag) )
		isoform = sapply( tag[range] , function(x) paste( sort(subset(LRcorres,correctTag%in%x)$SRindex) , collapse="_" ) )
		
		correctIsoformExp = tapply( read$coverage[range] , isoform , sum )
		srindex = lapply( strsplit( names(correctIsoformExp) , "_" ) , function(x) as.integer(x) )
		correctIsoform = lapply( srindex , function(ind) { 
			x=src[ ind , c('chr','start','end') ]
			x[ order(x$start), ]
			} )
		
		
		lrCount = tapply( LRcorres$count , LRcorres$SRindex , sum )
		LRjuncCount = data.frame( lrCount , src[as.integer(names(lrCount)),] )
		colnames(LRjuncCount) = c("lrCount","srCount","chr","start","end","motif")
		
		correctIsoformNum = length(unique(isoform))
		correctIsoformFrac = sum(correctIsoformExp)/sum(read$coverage) 
		
		list( isoform=correctIsoform , num=correctIsoformNum , exp=as.numeric(correctIsoformExp) , frac=correctIsoformFrac , LRjuncCount=LRjuncCount )
		
	}

	names(allCorrectIsoform) = CHR
  
	allCorrectIsoform
	
}
