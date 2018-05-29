echo "keyword: $1"
echo "Input Sam: $2"
echo "Output: $3"

samtools view -H $2 > $3
grep $1 $2 >> $3

