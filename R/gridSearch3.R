#' gridSearch
#'
#' @param LRjunc 
#' @param SRjunc 
#' @param thresSR 
#' @param thresDis 
#' @param adjustNCjunc 
#' @param matchedLS 
#'
#' @return
#' @export
#'
#' @examples
gridSearch = function( LRjunc , SRjunc , thresSR=c(1:30) , thresDis=c(1:20) , adjustNCjunc=TRUE , matchedLS=NULL )
{
  
  library(nnls)
  
  CHR = intersect( names(SRjunc) , names(LRjunc) )
  
  P = foreach( i = 1:length(thresSR) ) %dopar%
  {
    res = c()
    for( j in 1:length(thresDis) )
    {	
      ts = thresSR[i]
      td = thresDis[j]
      
      #cat(ts,td,"\n")
      
      LSjuncCount = foreach( k=1:length(CHR) ) %dopar%
      {
        #cat(k,'\n')
        chr = CHR[k]
        lrc = LRjunc[[chr]]
        src = SRjunc[[chr]]
        
        if( !is.null(matchedLS) )
        {
          SRmatch = matchedLS[[chr]]
          range = which( SRmatch[,'LSdis']<=td & src$count[SRmatch[,'SRindex']]>=ts ) 
          LRcorres = data.frame(lrc,SRmatch)[range, ]
          lrCount = tapply( LRcorres$count , LRcorres$SRindex , sum )
          LSjc = data.frame( lrCount , src[as.integer(names(lrCount)),] )
        } else {
          
          src = subset( src , count>=ts )
          if( nrow(src)>0 & nrow(lrc)>0 )
          {
            SRmatch = matchLSjuncOneChr( lrc , src )
            range = which( SRmatch[,'LSdis']<=td ) 
            LRcorres = data.frame(SRmatch,lrc)[range, ]
            lrCount = apply( as.matrix(LRcorres[,-c(1:5) ]) , 2 , 
                             function(x) tapply( x , LRcorres$SRindex , sum ) )
            LSjc = data.frame( src[as.integer(rownames(lrCount)),], lrCount )
          } else {
            LSjc = NULL
          }
          
        }
        
        LSjc
        
      }
      
      LSjuncCount = do.call(rbind,LSjuncCount)
      colnames(LSjuncCount)[1] = c("srCount")
      ind = which( rowSums(LSjuncCount[,-c(1:5)]) > 0)
      
      if(ncol(LSjuncCount)>6)
      {
        y = log(LSjuncCount[,1]+1)
        x = as.matrix( log(LSjuncCount[,-c(1:5)]+1) )
        
        # z = lm( y ~ x ) 
        # fit = predict(z)
        
        posCoef = as.matrix( coef(nnnpls(x,y,con=rep(1,ncol(x)))) )
        fit = x %*% as.matrix(posCoef)

      } else {
        y = LSjuncCount[,1]
        fit = log(LSjuncCount[,6]+1)
      }
      
      if( adjustNCjunc )
      {
        fracNCjunc = 1 - mean( LSjuncCount$motif[ind]==0 )
        res[j] = cor.test( fit , y , method="spearman" )$estimate * fracNCjunc
        
      } else {
        res[j] = cor.test( fit , y , method="spearman" )$estimate 
      }
      #res[j] = KL.Dirichlet( LSjuncCount[ind,1] , LSjuncCount[ind,2], a1=1/length(ind), a2=1/length(ind) )
      #res[j] = mi.Dirichlet( t(LSjuncCount[ind,1:2]) , a=1/(2*length(ind) )  )
      #res[j] = cor.test( log2(LSjuncCount[ind,1]+1) , log2(LSjuncCount[ind,2]+1) , method="pearson" )$estimate
      
      cat("ts=",ts,"td=",td,"score=",res[j],"\n")
      
    }
    res
  }
  
  P = do.call(rbind,P)
  P
}

