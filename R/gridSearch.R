#' gridSearch
#'
#' @param LRjunc 
#' @param SRjunc 
#' @param thresSR 
#' @param thresDis 
#' @param matchedLS 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
gridSearch = function( LRjunc , SRjunc , thresSR=c(1:30) , thresDis=c(1:20) , matchedLS=NULL )
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
      
      LRjuncCount = foreach( k=1:length(CHR) ) %dopar%
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
      LRjuncCount = do.call(rbind,LRjuncCount)
      colnames(LRjuncCount) = c("lrCount","srCount","chr","start","end","motif")
      ind = which(LRjuncCount$lrCount>0)
      res[j] = cor.test( LRjuncCount[ind,1] , LRjuncCount[ind,2] , method="spearman" )$estimate
      #res[j] = KL.Dirichlet( LRjuncCount[ind,1] , LRjuncCount[ind,2], a1=1/length(ind), a2=1/length(ind) )
      #res[j] = mi.Dirichlet( t(LRjuncCount[ind,1:2]) , a=1/(2*length(ind) )  )
      #res[j] = cor.test( log2(LRjuncCount[ind,1]+1) , log2(LRjuncCount[ind,2]+1) , method="pearson" )$estimate
    }
    res
  }
  
  P = do.call(rbind,P)
  ij = which( P==max(P) , arr.ind=T )
  ij
  
  P
}

