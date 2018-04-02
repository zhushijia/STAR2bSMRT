library(polyester)
library(Biostrings)
library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
library(foreach)
library(doMC)


outputDir = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/shijia/"
fasta_file = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/641/Exp/isoform_ts5_td17.fa"
readlen = 50


fasta = readDNAStringSet(fasta_file)
exp = sapply( strsplit( names(fasta) , "_" ) , function(x) as.integer(gsub("exp","",x[2])) )
fold_changes = matrix( rep(1,2*length(fasta)) , nrow=length(fasta) )
head(fold_changes)

# subset the FASTA file to first 20 transcripts
setwd(outputDir)
writeXStringSet( fasta, 'Isoseq_isoforms.fa')
writeXStringSet( fasta[exp>quantile(exp,0.5)], 'Isoseq_isoforms_highExp.fa')
writeXStringSet( fasta[exp<quantile(exp,0.5)], 'Isoseq_isoforms_lowExp.fa')

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round( exp * width(fasta) / readlen )

# simulation call:
simulate_experiment('Isoseq_isoforms.fa', readlen=readlen , reads_per_transcript=readspertx, 
                    num_reps=c(1,1), fold_changes=fold_changes, outdir='simulated_reads') 





module load kallisto
cd /sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/shijia/simulated_reads
perl ../fasta_to_fastq.pl sample_01_1.fasta > sample_01_1.fastq
perl ../fasta_to_fastq.pl sample_01_2.fasta > sample_01_2.fastq

fasta="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/shijia/Isoseq_isoforms.fa"
kallisto index -i transcripts.idx $fasta
#kallisto quant -i transcripts.idx -o sample_01 -b 10 --pseudobam sample_01_1.fastq sample_01_2.fastq > pseudo.sam
kallisto quant -i transcripts.idx -o sample_01 -b 10 sample_01_1.fastq sample_01_2.fastq 

fasta="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/shijia/Isoseq_isoforms_highExp.fa"
kallisto index -i transcripts.idx $fasta
kallisto quant -i transcripts.idx -o sample_01_high -b 10 sample_01_1.fastq sample_01_2.fastq

fasta="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/shijia/Isoseq_isoforms_lowExp.fa"
kallisto index -i transcripts.idx $fasta
kallisto quant -i transcripts.idx -o sample_01_low -b 10 sample_01_1.fastq sample_01_2.fastq




setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/shijia/simulated_reads/sample_01_high")
x=read.table("abundance.tsv",header=T)
exp = sapply( strsplit( as.character(x[,1]) , "_" ) , function(x) as.integer(gsub("exp","",x[2])) )
tpm = x[,5]
cor(exp,tpm)
cor(log(exp),log(tpm))




alignments="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/641/SR/alignments.bam"
gff="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/641/Exp/isoform_ts5_td17.gff"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/641/SR/"
. /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/STAR2bSMRT/data/stringTie.sh $alignments $gff $outputDir

