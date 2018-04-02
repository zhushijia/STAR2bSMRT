echo "stringTie isoform assembly ... ... "
echo "alignment bam or sam: $1"
echo "annotation gtf or gff: $2"
echo "outputDir: $3"

if [ ! -d $3 ]; then 
mkdir -p $3
fi
cd $3

module load stringtie
module load samtools
cd $outputDir
samtools view -Su $1 | samtools sort - $1.sorted
stringtie $1.sorted.bam -G $2 -o stringtie.gtf -p 4 -v -C stringtie.coverage -A stringtie.abundance

