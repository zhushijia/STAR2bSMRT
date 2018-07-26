setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/minimap2/SRR5286959_1")



bed = read.table("SRR5286959.sort.bed12.bed",sep="\t",header=T)
colnames(bed) = c("chrom","start","end","name","score","strand",
"thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")

f = function(y) sapply( strsplit( y , "," ) , function(x) as.integer(x) )
blockSizes = lapply( as.character(bed$blockSizes) , f )
blockStarts = lapply( as.character(bed$blockStarts) , f )

junc = sapply( 1:nrow(bed) , function(i) {
	exonStart = bed$start[i] + blockStarts[[i]]
	exonEnd = exonStart + blockSizes[[i]]
	paste( paste( exonEnd[-length(exonEnd)]+1 , exonStart[-1] , sep=",") , collapse="," )
} )

read = data.frame(id=bed$name, chr=bed$chrom , strand=bed$strand , start=bed$start, end=bed$end , junc=junc )
write.table(read,'alignments.read',sep='\t',quote=F,col.names=TRUE,row.names=FALSE)




1485,6337,6474,6560,6814,7552,7815,10282,10367,10647
1485,6337,6474,6560,6814,7552,7815,10282,10367,10647