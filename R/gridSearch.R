
gridSearch = function( LRjunc , SRjunc , matchedLS , thresSR=c(1:30) , thresDis=c(1:20) )
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
				SRmatch = matchedLS[[chr]]
				
				range = which( SRmatch[,'LSdis']<=td & src$count[SRmatch[,'SRindex']]>=ts ) 
				LRcorres = data.frame(lrc,SRmatch)[range, ]

				lrCount = tapply( LRcorres$count , LRcorres$SRindex , sum )
				data.frame( lrCount , src[as.integer(names(lrCount)),] )
			}
			LRjuncCount = do.call(rbind,LRjuncCount)
			colnames(LRjuncCount) = c("lrCount","srCount","chr","start","end","motif")
			ind = which(LRjuncCount$lrCount>0)
			res[j] = cor.test( LRjuncCount[ind,1] , LRjuncCount[ind,2] , method="spearman" )$estimate
		}
		res
	}
	
	P = do.call(rbind,P)
	ij = which( P==min(P) , arr.ind=T )
	ij
	
	P
}

