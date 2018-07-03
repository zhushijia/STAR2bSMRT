annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1ExonAnnotations.txt',sep='\t',header=T)
annotation[2,2] = annotation[3,3]+1
annotation[25,3] = annotation[24,2]-1

######################################
annotation[27,] = annotation[8,]
annotation[8,3] = 50848363
annotation[27,2] = 50848364
annotation[,1] = as.character(annotation[,1])
annotation[27,1] = "exon7c"
annotation[27,4] = 23


######################################
annotation[27,1] = "exon7a"
annotation[7,1] = "exon15_16"
annotation = annotation[c(1:6,27,8:16,17:26),]

annotation$fullStop = annotation$Stop
annotation$fullStart = annotation$Start
rownames(annotation) = 1:nrow(annotation)

annotation$fullStop[2] = annotation$fullStop[3]
annotation$fullStart[8] = annotation$fullStart[7]
annotation$fullStart[25] = annotation$fullStart[24]

len = ( annotation$fullStart - annotation$fullStop )

write.table(annotation, '/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep="\t",col.names=T,row.names=F,quote=F)


