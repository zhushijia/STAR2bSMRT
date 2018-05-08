system( "samtools view -bS Aligned.out.sam | bamToBed -i > Aligned.out.bed" )
library(data.table)
library(Biostrings)

bed = fread("Aligned.out.bed")
FL = bed[ V1=='chr17' & V2<90037384 & V3>91087954 ]

setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN")
fa = readDNAStringSet("SRR1184043_4.fa")
FLfa = fa[ names(fa) %in% as.character(FL$V4) ]
writeXStringSet(FLfa,"SRR1184043_4_NRXN1_FL_chr17_90037384_91087954.fa")
