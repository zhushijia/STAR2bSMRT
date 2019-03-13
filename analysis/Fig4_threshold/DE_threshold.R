setwd("/Users/shijiazhu/Dropbox/Projects/NRXN/results/Figs")
de = read.table("DEGseq.txt",sep="\t",header=T)
de2 = subset(de,flncDEL1>=7 | flncDEL2>=7 | flncCont1>=7 | flncCont2>=7)
pairs(log(de2[,1:4]+1),pch=16,col='red')
pairs( de2[,1:4],pch=16,col='red')
write.table(de2,"DEGseq_threshold7.txt",sep="\t",col.names=T,row.names=F,quote=F)

sum( de2$fdr<0.01 & de2$log2FoldChange>log2(2)  )
sum( de2$fdr<0.01 & de2$log2FoldChange>log2(4)  )
sum( de2$fdr<0.01 & de2$log2FoldChange>log2(10)  )

sum( de2$fdr<0.01 & de2$log2FoldChange<(-log2(2))  )
sum( de2$fdr<0.01 & de2$log2FoldChange<(-log2(4))  )
sum( de2$fdr<0.01 & de2$log2FoldChange<(-log2(10))  )


X = de2[,1:4]
X[X<7] = 0
pairs( log10(X+1) ,pch=16,col='red')

