echo "Input: $1"
echo "Output: $2"

if echo $1 | grep "sam"; then
  samtools view -bS $1 | bedtools bamtobed -bed12 -i > $2.bed12
fi

if echo $1 | grep "bam"; then
  bedtools bamtobed -bed12 -i $1 > $2.bed12
fi

echo -e "count\tchr\tstart\tend" > $2

cat $2.bed12 | awk '{

if ( $10 > 1 )
{
  split($11,size,",")
  split($12,start,",")

  js = $2 + start[1] + size[1] + 1
  je = $2 + start[2]
  print $1 "\t" js "\t" je
  
  for ( i=3; i<=$10; i++ )
  {
    js = $2 + start[i-1] + size[i-1] + 1
    je = $2 + start[i]
    print $1 "\t" js "\t" je
  }
  
}
  
}'  | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' >> $2

