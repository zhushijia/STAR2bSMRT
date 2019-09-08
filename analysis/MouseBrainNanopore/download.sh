PRJEB25574

#Illumina HiSeq 4000 paired end sequencing; Mus musculus_Illumina_RB
module load sratoolkit/2.8.0
cd /hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MouseBrainNanopore
fastq-dump --split-files ERR2680425
gzip ERR2680425_1.fastq &
gzip ERR2680425_2.fastq &

#Nanopore RNAseq mouse brain
module load sratoolkit/2.8.0
cd /hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/MouseBrainNanopore
fastq-dump --split-files ERR2401483 &


wget http://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-DirectRNA.pass.dedup.NoU.fastq.hg38.minimap2.sorted.bam &
wget http://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-DirectRNA.pass.dedup.NoU.fastq.hg38.minimap2.sorted.bam.bai
wget http://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-DirectRNA.fail.dedup.NoU.fastq.hg38.minimap2.sorted.bam
wget http://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-DirectRNA.fail.dedup.NoU.fastq.hg38.minimap2.sorted.bam.bai


