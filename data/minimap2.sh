echo "minimap2 mapping ... ... "
echo "genomeDir: $1"
echo "LR: $2"
echo "outputDir: $3"

if [ ! -d $3 ]; then 
mkdir -p $3
fi
cd $3


module load minimap2/r572

#minimap2 -d ref.mmi ref.fa                     # indexing
#minimap2 -a ref.mmi reads.fq > alignment.sam   # alignment
minimap2 -ax map-pb  ref.fa pacbio-reads.fq > aln.sam   # for PacBio subreads
minimap2 -ax map-ont ref.fa ont-reads.fq > aln.sam      # for Oxford Nanopore reads
minimap2 -ax splice -uf -C5 ref.fa iso-seq.fq > aln.sam      # PacBio Iso-seq/traditional cDNA
minimap2 -ax splice ref.fa nanopore-cdna.fa > aln.sam        # Nanopore 2D cDNA-seq
minimap2 -ax splice -uf -k14 ref.fa direct-rna.fq > aln.sam  # Nanopore Direct RNA-seq
minimap2 -ax splice --splice-flank=no SIRV.fa SIRV-seq.fa    # mapping against SIRV control



LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/intermediate_files/isoseq_nfl.fastq"
ref="/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/gencode/GRCh37.liftover.from.GRCh38.p10.Release26/genome/GRCh37.primary_assembly.genome.fa"
outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Alzheimer_IsoSeq_2016/intermediate_files/nfl_minimap2"
cd $outputDir
#minimap2 -N0 -t10 -ax map-pb $ref $LR > aln.sam
/hpc/users/zhus02/schzrnas/sjzhu/bitbucket/minimap2/minimap2 --secondary=no -t 30 -ax splice -uf -C5 $ref $LR > aln.sam      # PacBio Iso-seq/traditional cDNA


samtools view -hF 4 aln.sam -o Align.out.sam
samtools view -b -hf 4 aln.sam -o unmapped.bam
samtools sort -n unmapped.bam unmapped.sort
bedtools bamtofastq -i unmapped.sort.bam -fq Unmapped.out.mate1



