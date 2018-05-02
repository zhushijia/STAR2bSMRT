calculateExonExp = function(cc_fracs,ccExp)
{
  exonExpExpEachTranscript = apply( cc_fracs , 2 , function(x) x*rowSums(ccExp) )
  colSums(exonExpExpEachTranscript)
}

calculateJuncExp = function(cc_fracs,ccExp)
{
  isoform = apply(cc_fracs , 1 , function(x) {
    iso1 = paste( which(x>0), collapse = "_" )
    paste0( "_" , iso1 , "_")
  })
  
  
  junc = matrix(ncol=ncol(cc_fracs),nrow=ncol(cc_fracs),data=-254)
  colnames(junc) = colnames(cc_fracs)
  rownames(junc) = colnames(cc_fracs)
  
  for(i in 1:ncol(cc_fracs))
  {
    for(j in i:ncol(cc_fracs))
    {
      if(i==j)
      {
        junc[i,i] = sum( grepl( paste("",i,"",sep="_") , isoform ) )
      } else {
        junc[i,j] = sum( grepl( paste("",i,j,sep="_") , isoform ) )
      }
    }
    
  }
  
 junc
}

calculateJuncFrac = function(cc_fracs,ccExp)
{
  isoform = apply(cc_fracs , 1 , function(x) {
    iso1 = paste( which(x>0), collapse = "_" )
    paste0( "_" , iso1 , "_")
  })
  
  
  junc = matrix(ncol=ncol(cc_fracs),nrow=ncol(cc_fracs),data=0)
  colnames(junc) = colnames(cc_fracs)
  rownames(junc) = colnames(cc_fracs)
  
  for(i in 1:ncol(cc_fracs))
  {
    for(j in i:ncol(cc_fracs))
    {
      if(i==j)
      {
        junc[i,i] = sum( grepl( paste("",i,"",sep="_") , isoform ) )
      } else {
        junc[i,j] = sum( grepl( paste("",i,j,sep="_") , isoform ) )
      }
    }
    
  }
  
  frac = junc/diag(junc)
  frac[lower.tri(frac)] = -1
  diag(frac) = 0
  frac
}



caseExonExp = calculateExonExp(case_fracs,caseExp)
contExonExp = calculateExonExp(cont_fracs,contExp)

caseJuncExp = calculateJuncExp(case_fracs,caseExp)
contJuncExp = calculateJuncExp(cont_fracs,contExp)

caseJuncFrac = calculateJuncFrac(case_fracs,caseExp)
contJuncFrac = calculateJuncFrac(cont_fracs,contExp)


pdf("case.exon_junc_exp.pdf")
barplot(caseExonExp,border=F,las=3)
color = list( redgreen(256) , topo.colors(256) , cm.colors(256) , terrain.colors(256) , rainbow(256), heat.colors(256)  )
heatmap.2(caseJuncExp, col=color[[3]], scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
heatmap.2(caseJuncFrac, col=color[[3]], scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
dev.off()


pdf("cont.exon_junc_exp.pdf")
barplot(contExonExp,border=F,las=3)
color = list( redgreen(256) , topo.colors(256) , cm.colors(256) , terrain.colors(256) , rainbow(256), heat.colors(256)  )
heatmap.2(contJuncExp, col=color[[3]], scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
heatmap.2(contJuncFrac, col=color[[3]], scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
dev.off()


