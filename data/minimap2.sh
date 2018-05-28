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
