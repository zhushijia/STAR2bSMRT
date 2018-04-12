
module load minimap2/r572
ref="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a/SIRV_isoforms_multi-fasta_170612a.fasta"
seq="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/E2_Nanopore/SRR5286959.fastq"
cd /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/SRR5286959_1
minimap2 -ax splice --splice-flank=no $ref $seq > SRR5286959.sam
minimap2 -x ava-ont $seq $seq > SRR5286959.paf



module load star  
cd /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a
STAR --runMode genomeGenerate --genomeDir /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a/STAR_index --genomeFastaFiles SIRV_isoforms_multi-fasta_170612a.fasta  --runThreadN 10


cd /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/starLong/SRR5286959

seq="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/E2_Nanopore/SRR5286959.fastq"
genomeDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a/STAR_index"

/hpc/users/zhus02/schzrnas/sjzhu/bitbucket/STAR/STAR/source/STARlong \
--runMode alignReads \
--readNameSeparator space \
--outFilterMultimapScoreRange 1 \
--outFilterMismatchNmax 2000 \
--scoreGapNoncan -20 \
--scoreGapGCAG -4 \
--scoreGapATAC -8 \
--scoreDelOpen -1 \
--scoreDelBase -1 \
--scoreInsOpen -1 \
--scoreInsBase -1 \
--alignEndsType Local \
--seedSearchStartLmax 50 \
--seedPerReadNmax 100000 \
--seedPerWindowNmax 1000 \
--alignTranscriptsPerReadNmax 100000 \
--alignTranscriptsPerWindowNmax 10000 \
--outReadsUnmapped Fastx \
--outFilterIntronStrands None \
--outSAMattributes NH HI AS nM NM MD jM jI \
--genomeDir $genomeDir \
--readFilesIn $seq


cd /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a
library(data.table)
file = "SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf"
gtf = fread(file)
gtf$V9 = sapply( strsplit(gtf[,V9],"; ") , function(x) paste0(x[1],"; ",x[2],";") )

transcript = sapply( strsplit(gtf[,V9]," ") , function(x) gsub( '"|;' , "" , x[4]) )
exon1 = split(gtf,transcript)

exon2 = list()
for(i in 1:length(exon1))
{
	exon = exon1[[i]]
	exon = exon[order(exon[,4]),]
	trans = data.table( exon[1,1:2] , "transcript" , exon[1,4] , exon[nrow(exon),5] , exon[1,6:9] )
	colnames(trans) = paste("V",1:9,sep="")
	exon2[[i]] = rbind(trans,exon)
}
exon3 = do.call(rbind,exon2)
write.table( exon3 , "SIRV_isoforms_multi-fasta-annotation_C_170612a.gff" , col.names=F , row.names=F , sep="\t" , quote=F )

library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
library(Biostrings)
setwd("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a")
genome = readDNAStringSet("SIRV_isoforms_multi-fasta_170612a.fasta")
gff = readGff("SIRV_isoforms_multi-fasta-annotation_C_170612a.gff",chrom=NULL,s=0,e=Inf)
seq = generateSeq( genome , isoform=gff )
writeXStringSet( seq$dna , "SIRV_isoforms_multi-fasta-annotation_C_170612a_shijia.fa" )
	























