annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1.txt',sep='\t')
setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping")
samples = dir()[!grepl("_0",dir())]
exon = lapply( samples,function(sample) read.table( paste0(sample,"/NRXN1.exon.txt"),sep="\t",header=F) ) 
range = match( paste(annotation$V2,annotation$V3) , paste(exon[[1]]$V3,exon[[1]]$V4) )
exp = lapply( exon,function(x) x[range,])
names(exp) = samples

x = do.call(cbind,lapply(exp,function(x)x[,7]))
rownames(x) = c(24:1)
info = data.frame(exp[[1]][,1:6],x)
write.table(info,"NRNX1_24ExonExp.txt",sep="\t",col.names=T,row.names=F,quote=F)

require('gplots')
y = apply(x,2,function(z) z/mean(z) )
pdf("NRNX1_24ExonExp_heatmap.pdf")
heatmap.2(y, col=heat.colors(256), scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
dev.off()



