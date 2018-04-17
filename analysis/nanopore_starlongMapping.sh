module load star  
cd /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a
STAR --runMode genomeGenerate --genomeDir /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a/STAR_index --genomeFastaFiles SIRV_isoforms_multi-fasta_170612a.fasta  --runThreadN 10
STAR --runMode genomeGenerate --genomeDir /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a/STAR_index2 --genomeFastaFiles SIRV_isoforms_multi-fasta_170612a.fasta  --runThreadN 10 --genomeSAindexNbases 8


samples="SRR5286956_Lexogen_100fg_2_Smartseq2_ONT SRR5286957_Lexogen_100fg_1_Smartseq2_ONT SRR5286958_Lexogen_10fg_2_Smartseq2_ONT"


for sample in $samples; do

echo $sample

outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/starLong/${sample}"
seq="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/E2_Nanopore/${sample}.fastq"
genomeDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a/STAR_index"

mkdir $outputDir
cd $outputDir

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

done




sample="SRR6058583_NatMethod_E2"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/starLong/${sample}"
seq="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/E2_Nanopore/${sample}.fastq"
genomeDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a/STAR_index2"

mkdir $outputDir
cd $outputDir

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




module load minimap2/r572
ref="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/reference/SIRV_Set1_Sequences_170612a/SIRV_isoforms_multi-fasta_170612a.fasta"
seq="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/E2_Nanopore/SRR6058583_NatMethod_E2.fastq"
cd /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/SIRV/result/E2_Nanopore/minimap2/SRR6058583_NatMethod_E2
minimap2 -ax splice --splice-flank=no $ref $seq > SRR6058583_NatMethod_E2.sam






















