parent="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping"
setwd(parent)
samples=dir()[ grepl("^EF",dir()) ]

library(foreach)
library(doMC)
registerDoMC(10)

system.time( res <- foreach( i = 1:length(samples) ) %dopar%
{
	sample = samples[i]
	bam = paste0(parent,"/",sample,"/alignments.bam")
	bedFile = paste0(parent,"/NRXN1_hg19_ExonAnnotations_shijia.bed")
	outputFile = paste0(parent,"/",sample,"/NRXN1_ExonExp_bedTools.txt")
	cmd = paste0( "bedtools multicov -bams " , bam , " -bed " , bedFile , " > " , outputFile )
	message(cmd)
	system(cmd)
} )

#########################

parent="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping"
setwd(parent)

getCount = function( ef )
{
  bed = read.table( paste0(parent,"/",ef,"/NRXN1_ExonExp_bedTools.txt"),sep='\t',header=F)
  count = bed$V5
  names(count) = as.character(bed$V4)
  len = bed$V3 - bed$V2 + 1
  count = count / len
  count['exon3b'] = count['exon3b'] - count['exon3a']
  count['exon7a'] = count['exon7a'] - count['exon7b']
  count['exon23a'] = count['exon23a'] - count['exon23b']
  
  count['exon3b'] = count['exon3b'] - count['exon3a']
  count['exon7a'] = count['exon7a'] - count['exon7b']
  count['exon23a'] = count['exon23a'] - count['exon23b']
  
  
  geneCount = sum(count) / sum(len)
  
  exonCount/geneCount
  
  
}

getPSI = function( count )
{
  psi = count/mean(count)
  SS = c('exon3a','exon3b','exon4','exon5','exon7a','exon7b','exon12','exon17','exon21','exon23a','exon23b')
  psi[ names(psi) %in% SS ]
}


PBS_2607 = getCount('EF15')
KCL_2607 = getCount('EF17')
PBS_553 = getCount('EF21')
KCL_553 = getCount('EF23')
PBS_641 = getCount('EF18')
KCL_641 = getCount('EF20')

KCL_case = KCL_641
PBS_case = PBS_641
KCL_cont = KCL_2607 + KCL_553
PBS_cont = PBS_2607 + PBS_553

ratio_case = getPSI(KCL_case)/getPSI(PBS_case)
ratio_cont = getPSI(KCL_cont)/getPSI(PBS_cont)

PSI_ratio = c( matrix( rbind(ratio_case,ratio_cont,0) , nrow=1 ))


barplot(PSI_ratio,col=rep(c('darkred','grey','black'),11),las=3)
dev.off()


