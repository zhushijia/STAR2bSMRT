library(Biostrings)
setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN")
fa1 = readDNAStringSet("SRR1184043.fa")
fa2 = readDNAStringSet("SRR1184044.fa")
names(fa1) = paste("SRR1184043",names(fa1),sep="_")
names(fa2) = paste("SRR1184044",names(fa2),sep="_")
fa12 = c(fa1,fa2)

writeXStringSet( fa12 , "SRR1184043_4.fa" )
