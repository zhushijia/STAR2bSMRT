echo "1pass STARlong mapping ... ... "
echo "genomeDir: $1"
echo "LR: $2"
echo "outputDir: $3"
echo "cores: $4"


if [ ! -d $3 ]; then 
mkdir -p $3
fi
cd $3

/sc/orga/projects/schzrnas/sjzhu/bitbucket/STAR/STAR/source/STARlong \
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
--runThreadN $4 \
--genomeDir $1 \
--readFilesIn $2


# genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
# LR="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Sample_553-x-C-FBN/myoutputdir/corrected_LR.fa"
# outputDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Sample_553-x-C-FBN/starLong/corrected"

