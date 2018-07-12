parent="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping"
setwd(parent)
samples=dir()[ grepl("^EF",dir()) ]

library(foreach)
library(doMC)
registerDoMC(20)

system.time( res <- foreach( i = 1:length(samples) ) %dopar%
{
	sample = samples[i]
	bam = paste0(parent,"/",sample,"/alignments.bam")
	bedFile = paste0(parent,"/NRXN1_hg19_fullExonAnnotations_shijia.bed")
	outputFile = paste0(parent,"/",sample,"/NRXN1_fullExonExp_bedTools.txt")
	cmd = paste0( "bedtools multicov -bams " , bam , " -bed " , bedFile , " > " , outputFile )
	message(cmd)
	system(cmd)
} )

#########################

parent="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping"
setwd(parent)

getCount = function( ef )
{
  bed = read.table( paste0(parent,"/",ef,"/NRXN1_fullExonExp_bedTools.txt"),sep='\t',header=F)
  count = bed$V5
  len = bed$V3 - bed$V2 + 1
  names(count) = as.character(bed$V4)
  names(len) = as.character(bed$V4)
  
  count['exon3a'] = count['exon3a'] - count['exon3b']
  count['exon7b'] = count['exon7b'] - count['exon7a']
  count['exon23b'] = count['exon23b'] - count['exon23a']
  
  len['exon3a'] = len['exon3a'] - len['exon3b']  + 1
  len['exon7b'] = len['exon7b'] - len['exon7a'] + 1
  len['exon23b'] = len['exon23b'] - len['exon23a'] + 1

  exonCount = count / len
  geneCount = sum(count[-c(1,26)]) / sum(len[-c(1,26)])
  
  psi = exonCount/geneCount
  SS = c('exon3a','exon3b','exon4','exon5','exon7a','exon7b','exon12','exon17','exon21','exon23a','exon23b')
  psi[ names(psi) %in% SS ]
}

######################################################
PBS_2607 = getCount('EF15')
KCL_2607 = getCount('EF17')
PBS_553 = getCount('EF21')
KCL_553 = getCount('EF23')
PBS_641 = getCount('EF18')
KCL_641 = getCount('EF20')

KCL_case = KCL_641
PBS_case = PBS_641
KCL_cont = ( KCL_2607 + KCL_553 ) / 2
PBS_cont = ( PBS_2607 + PBS_553 ) / 2

ratio_case = KCL_case/PBS_case
ratio_cont = KCL_cont/PBS_cont

PSI_ratio = c( matrix( rbind(ratio_cont,ratio_case,0) , nrow=1 ))
names(PSI_ratio) = c( matrix(rbind(names(ratio_case),"",""), nrow=1 ))

pdf("PSI_FC_KCLvsPSB.pdf")
barplot(PSI_ratio,col=rep(c('grey','darkred','black'),11),las=3)
dev.off()

######################################################
PBS_Neuron_w2 = ( getCount('EF27') + getCount('EF31') ) / 2
PBS_Neuron_w4 = ( getCount('EF28') + getCount('EF32') ) / 2
PBS_Neuron_w6 = ( getCount('EF15') + getCount('EF21') ) / 2
PBS_NPC = getCount('EF25') 

ratio_w2 = PBS_Neuron_w2/PBS_NPC
ratio_w4 = PBS_Neuron_w4/PBS_NPC
ratio_w6 = PBS_Neuron_w6/PBS_NPC

PSI_ratio = c( matrix( rbind(ratio_w2,ratio_w4,ratio_w6,0) , nrow=1 ))
names(PSI_ratio) = c( matrix(rbind(names(ratio_w2),"","",""), nrow=1 ))

pdf("PSI_FC_NeuronWeek2_4_6vsNPC.pdf")
barplot(PSI_ratio,col=rep(c('grey','darkblue','darkred','black'),11),las=3)
dev.off()









