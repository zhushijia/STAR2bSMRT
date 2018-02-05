echo "transcripts: $1"
echo "R1: $2"
echo "R2: $3"
echo "OutputDir: $4"

if [ ! -d $4 ]; then 
mkdir -p $4
fi
cd $4

salmon index -t $1 -i transcripts_index
salmon quant -i transcripts_index -l A -1 $2 -2 $3 -p 8 -o quants
