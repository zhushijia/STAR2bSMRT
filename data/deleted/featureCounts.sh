echo "calling read counts ... ... "
echo "gtf: $1"
echo "bam: $2"
echo "outputDir: $3"
echo "feature: $4"
echo "cores: $5"

module load subread

if [ ! -d $3 ]; then 
mkdir -p $3
fi
cd $3

if [ $4 == 'exon' ]; then
echo "... generating exonReadCount.txt"
featureCounts -T $5 -t exon -f -O -p -a  $1  -o  exonReadCount.txt  $2
fi

if [ $4 == 'gene' ]; then
echo "... generating primaryGeneReadCount.txt"
featureCounts -T $5 -t exon -g gene_id --primary -O -p -a $1 -o primaryGeneReadCount.txt $2
fi

if [ $4 == 'transcript' ]; then
echo "... generating transcriptReadCount.txt"
featureCounts -T $5 -t exon -g transcript_id -p -a $1 -o transcriptReadCount.txt $2
fi


#if [ $4 == 'exon' ]; then
#featureCounts -T $5 -t exon -g gene_id -p -a $1 -o geneReadCount.txt $2
#fi
