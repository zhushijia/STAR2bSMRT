#' generateCorrectedIsoform
#'
#' @param LRjunc 
#' @param SRjunc 
#' @param LRtag 
#' @param LRread 
#' @param ts 
#' @param td 
#' @param matchedLS 
#' @param fuzzyMatch
#'
#' @return
#' @export
#'
#' @examples
generateCorrectedIsoform = function( LRjunc , SRjunc, LRtag , LRread  , ts , td , matchedLS=NULL , fuzzyMatch=100 )
{
  
  CHR = intersect( names(SRjunc) , names(LRjunc) )
  
  allCorrectIsoform = foreach( k=1:length(CHR) ) %dopar% {
    
    chr = CHR[k]
    lrc = LRjunc[[chr]]
    src = SRjunc[[chr]]
    tag = LRtag[[chr]]
    read = LRread[[chr]]
    
    if( !is.null(matchedLS) )
    {
      SRmatch = matchedLS[[chr]]
      index = which( SRmatch[,'LSdis']<=td & src$count[SRmatch[,'SRindex']]>=ts ) 
    } else {
      src = subset( src , count>=ts )
      SRmatch = matchLSjuncOneChr( lrc , src ,fuzzyMatch )
      index = which( SRmatch[,'LSdis']<=td ) 
    }
    
    #LRcorres = data.frame(lrc[index,],SRmatch[index,])
    LRcorres = data.frame(SRmatch,lrc)[index, ]
    correctTag = paste(LRcorres$start , LRcorres$end , sep="," )
    range = sapply( tag , function(x) all(x%in%correctTag) )
    isoform = sapply( tag[range] , function(x) paste( sort(subset(LRcorres,correctTag%in%x)$SRindex) , collapse="_" ) )
    
    #correctIsoformExp = tapply( read$coverage[range] , isoform , sum )
    #srindex = lapply( strsplit( names(correctIsoformExp) , "_" ) , function(x) as.integer(x) )
    
    info = data.frame( read[range,] , isoform  )
    correctIsoformExp = cast( info , isoform ~ group, sum , value="coverage")
    srindex = lapply( strsplit( as.character(correctIsoformExp$isoform) , "_" ) , 
                      function(x) as.integer(x) )
    
    correctIsoform = lapply( srindex , function(ind) { 
      x=src[ ind , c('chr','start','end') ]
      x[ order(x$start), ]
    } )
    
    
    #lrCount = tapply( LRcorres$count , LRcorres$SRindex , sum )
    #LSjuncCount = data.frame( lrCount , src[as.integer(names(lrCount)),] )
    #colnames(LSjuncCount) = c("lrCount","srCount","chr","start","end","motif")
    
    lrCount = apply( as.matrix(LRcorres[,-c(1:5) ]) , 2 , 
                     function(x) tapply( x , LRcorres$SRindex , sum ) )
    LSjuncCount = data.frame( src[as.integer(rownames(lrCount)),], lrCount )
    colnames(LSjuncCount)[1] = c("srCount")
    
    correctIsoformNum = length(unique(isoform))
    
    if( ncol(correctIsoformExp) == 2 )
    {
      correctIsoformFrac = sum(correctIsoformExp[,2])/sum(read$coverage)
    }
     
    list( isoform=correctIsoform , num=correctIsoformNum , exp=correctIsoformExp , frac=correctIsoformFrac , LSjuncCount=LSjuncCount )
    
  }
  
  names(allCorrectIsoform) = CHR
  
  allCorrectIsoform
  
}
