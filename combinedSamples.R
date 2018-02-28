library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/"
setwd(folder)
samples= dir()
samples = samples[samples!="combinedSamples"]

gffs = lapply( samples , function(sample)
{
  files = dir( paste0(sample,"/Exp") )
  gff = paste0(sample,"/Exp/",files[grep("gff",files)])
  readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(gffs) = samples


setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/combinedSamples")
allgff = gffs[[1]]
for(i in 2:5)
{
  allgff = unionGff(allgff,gffs[[i]])
  cat(length(allgff),"")
}
writeGff( allgff , file = "allsamples.gff" )

genomeFasta = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
genome = readDNAStringSet(genomeFasta)
seq = generateSeqFromExon( genome , isoform=allgff )
fa = seq$dna
names(fa) = paste(1:length(fa),names(fa))
writeXStringSet( fa[which(seq$translated)] , "translatedIsoform.fa" )
writeXStringSet( fa[sapply(seq$dna,nchar)<5000] , "allIsoforms.fa" )


fastaName = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/combinedSamples/allIsoforms.fa"
folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/combinedSamples/allIsoforms/"

kallisto = list()
sample = "581"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/581/581.R2.fastq"
outputDir=paste0(folder,sample)
kallisto[[sample]] = kallistoQuant( fastaName , SR1 , SR2 , outputDir )

sample = "641"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/641/641.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/641/641.R2.fastq"
outputDir=paste0(folder,sample)
kallisto[[sample]] = kallistoQuant( fastaName , SR1 , SR2 , outputDir )

sample = "2607"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/2607/2607.R2.fastq"
outputDir=paste0(folder,sample)
kallisto[[sample]] = kallistoQuant( fastaName , SR1 , SR2 , outputDir )

sample = "553"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/553/553.R1.fastq"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/MiSeq/KM1707142-R1-44416635-unzip/553/553.R2.fastq"
outputDir=paste0(folder,sample)
kallisto[[sample]] = kallistoQuant( fastaName , SR1 , SR2 , outputDir )

sample = "642"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/642-2-redo_S14_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Cleaned_targetShortRead/CG1709012_R1/642-2-redo_S14_R2_001.fastq.gz"
outputDir=paste0(folder,sample)
kallisto[[sample]] = kallistoQuant( fastaName , SR1 , SR2 , outputDir )


gffOfInterest = allgff[seq$translated]
tpm = sapply(kallisto,function(x)x$tpm)
sapply( 1:length(kallisto) , function(i) {
  exp = sapply( strsplit(names(gffs[[i]]),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
  range = matchGff( gffs[[i]] , gffOfInterest )
  exp = exp[!is.na(range)]
  range = range[!is.na(range)]
  r = cor.test( log(exp+1) , log(tpm[range,i]+1) )$estimate
  p = cor.test( log(exp+1) , log(tpm[range,i]+1) )$p.val
  c(r,p,length(range))

} )



