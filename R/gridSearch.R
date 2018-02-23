#' gridSearch
#'
#' @param LRjunc 
#' @param SRjunc 
#' @param thresSR 
#' @param thresDis 
#' @param matchedLS 
#' @param adjustNCjunc a boolean value indicating whether or not to control the non-carnonical junction site fraction
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
gridSearch = function( LRjunc , SRjunc , thresSR=c(1:30) , thresDis=c(1:20) , adjustNCjunc=TRUE , matchedLS=NULL )
{
  
  CHR = intersect( names(SRjunc) , names(LRjunc) )
  
  P = foreach( i = 1:length(thresSR) ) %dopar%
  {
    res = c()
    for( j in 1:length(thresDis) )
    {	
      ts = thresSR[i]
      td = thresDis[j]
      
      cat(ts,td,"\n")
      
      LSjuncCount = foreach( k=1:length(CHR) ) %dopar%
      {
        chr = CHR[k]
        lrc = LRjunc[[chr]]
        src = SRjunc[[chr]]
        
        if( !is.null(matchedLS) )
        {
          SRmatch = matchedLS[[chr]]
          range = which( SRmatch[,'LSdis']<=td & src$count[SRmatch[,'SRindex']]>=ts ) 
        } else {
          src = subset( src , count>=ts )
          SRmatch = matchLSjuncOneChr( lrc , src )
          range = which( SRmatch[,'LSdis']<=td ) 
        }
        
        LRcorres = data.frame(lrc,SRmatch)[range, ]
        lrCount = tapply( LRcorres$count , LRcorres$SRindex , sum )
        data.frame( lrCount , src[as.integer(names(lrCount)),] )
      }
      LSjuncCount = do.call(rbind,LSjuncCount)
      colnames(LSjuncCount) = c("lrCount","srCount","chr","start","end","motif")
      ind = which(LSjuncCount$lrCount>0)
      
      if( adjustNCjunc )
      {
        fracNCjunc = 1 - mean( LSjuncCount$motif[ind]==0 )
        res[j] = cor.test( LSjuncCount[ind,1] , LSjuncCount[ind,2] , method="spearman" )$estimate * fracNCjunc
        
      } else {
        res[j] = cor.test( LSjuncCount[ind,1] , LSjuncCount[ind,2] , method="spearman" )$estimate 
      }
      #res[j] = KL.Dirichlet( LSjuncCount[ind,1] , LSjuncCount[ind,2], a1=1/length(ind), a2=1/length(ind) )
      #res[j] = mi.Dirichlet( t(LSjuncCount[ind,1:2]) , a=1/(2*length(ind) )  )
      #res[j] = cor.test( log2(LSjuncCount[ind,1]+1) , log2(LSjuncCount[ind,2]+1) , method="pearson" )$estimate
    }
    res
  }
  
  P = do.call(rbind,P)
  P
}

