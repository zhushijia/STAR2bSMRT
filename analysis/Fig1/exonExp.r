gtf="/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/gencode/GRCh37.liftover.from.GRCh38.p10.Release26/annotation/gencode.v26lift37.annotation.gtf"
for sample in 642 adult_dlPFC1_10 adult_dlPFC1_13 fetal NRXN_adult_dlPFC1_12; do
cd /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/$sample/SR
echo $sample
pwd
featureCounts -T 10 -t exon -f -O -p -a $gtf -o exon.txt alignments.bam
grep ENSG00000179915 exon.txt > NRXN1.exon.txt
done


annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1.txt',sep='\t')
setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/")
samples = dir()[!grepl("combinedSamples",dir())]
exon = lapply( samples,function(sample) read.table( paste0(sample,"/SR/NRXN1.exon.txt"),sep="\t",header=F) ) 
range = match( paste(annotation$V2,annotation$V3) , paste(exon[[1]]$V3,exon[[1]]$V4) )
exp = lapply( exon,function(x) x[range,])
names(exp) = samples

x = do.call(cbind,lapply(exp,function(x)x[,7]))
info = data.frame(exp[[1]][,1:6],x)
setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/result/pipeline_nonAdjustNCjunc")
write.table(info,"NRNX1_24ExonExp.txt",sep="\t",col.names=T,row.names=F)













