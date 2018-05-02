case_fracs2 = case_fracs
cont_fracs2 = cont_fracs

range = matchGff(contGff,caseGff)
range = range[!is.na(range)]
case_fracs2[range, ] = -case_fracs2[range, ]

range = matchGff(caseGff,contGff)
range = range[!is.na(range)]
cont_fracs2[range, ] = -cont_fracs2[range, ]



plotHeatmap = function(cc_fracs,ccExp)
{
  require('gplots')
	ccExp = log(ccExp+1)
	x = cc_fracs
	#x = cbind( (cc_fracs) , 0 , 1-ccExp/max(ccExp) )
	#x = x[order(x[,2]),]
	#heatmap( x ,Rowv=NA,Colv=NA , scale="none" , col=redgreen(256) )
	heatmap.2(x, col=cm.colors(256), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
}

pdf("case.heatmap_translated_ConsiderOverlap.pdf")
plotHeatmap( case_fracs2[caseSeq$translated,] , caseExp[caseSeq$translated] )
dev.off()

pdf("cont.heatmap_translated_ConsiderOverlap.pdf")
plotHeatmap( cont_fracs2[contSeq$translated,] , contExp[contSeq$translated] )
dev.off()







color = list( redgreen(256) , topo.colors(256) , cm.colors(256) , terrain.colors(256) , rainbow(256), heat.colors(256)  )


plotExonCoexpress = function(cc_fracs, metrics=c("correlation","Euclidean"))
{
#colnames(cc_fracs)=1:ncol(cc_fracs)
if( metrics=="correlation")
{
	sds = apply(cc_fracs,2,sd)
	cc_fracs = cc_fracs[,sds>0]
	heatmap( 1-cor(cc_fracs) ,Rowv=NA,Colv=NA , scale="none")
}
if( metrics=="Euclidean")
	#heatmap( as.matrix(dist(t(cc_fracs))) ,Rowv=NA,Colv=NA , scale="none")
  heatmap.2(as.matrix(dist(t(cc_fracs))), col=heat.colors(256), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
}

pdf("case.exonCoexpress_translated.pdf")
plotExonCoexpress( case_fracs[caseSeq$translated,] , "Euclidean" )
dev.off()

pdf("cont.exonCoexpress_translated.pdf")
plotExonCoexpress( cont_fracs[contSeq$translated,] , "Euclidean" )
dev.off()


