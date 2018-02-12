echo "transcripts: $1"
echo "R1: $2"
echo "R2: $3"
echo "OutputDir: $4"

if [ ! -d $4 ]; then 
mkdir -p $4
fi
cd $4

module load kallisto
kallisto index -i transcripts.idx $1
kallisto quant -i transcripts.idx -o output -b 10 $2 $3
