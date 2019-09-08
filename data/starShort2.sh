echo "2pass STARshort mapping ... ... "
echo "genomeDir: $1"
echo "SR: $2 $3"
echo "outputDir: $4"
echo "SJ.tab: $5"

if [ ! -d $4 ]; then 
mkdir -p $4
fi
cd $4

if echo $2 | grep "gz$"  
then
STAR --genomeDir $1 \
--sjdbFileChrStartEnd $5 \
--readFilesIn $2 $3 \
--readFilesCommand zcat \
--outSAMattributes NH HI AS nM NM MD jM jI \
--outReadsUnmapped Fastx \
--outFilterIntronStrands None \
--outSAMstrandField intronMotif --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate \
--runThreadN 10 > alignments.bam
else
STAR --genomeDir $1 \
--sjdbFileChrStartEnd $5 \
--readFilesIn $2 $3 \
--outSAMattributes NH HI AS nM NM MD jM jI \
--outReadsUnmapped Fastx \
--outFilterIntronStrands None \
--outSAMstrandField intronMotif --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate \
--runThreadN 10 > alignments.bam
fi

samtools index alignments.bam



# --outFilterMultimapNmax 1  # reads mapped to maximum loci
# --chimSegmentMin 15
# --chimJunctionOverhangMin 15 



# genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
# R1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_2607-2-1-FBN/Sample_2607-2-1-FBN.R1.fastq"
# R2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_2607-2-1-FBN/Sample_2607-2-1-FBN.R2.fastq"
# outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_2607-2-1-FBN/starShortNew"
