echo "Input: $1"
echo "Output: $2"

echo -e "id\tmovieName\tZMWnumber\tsubread_start_end" > $2


if echo $1 | grep -P "sam$|bam$"  
then
samtools view $1 | \
awk '{ 
	split($1,a,"/")
	print $1 "\t" a[1] "\t" a[2] "\t" a[3]
}' >> $2
fi

if echo $1 | grep -P "fa$|fq$|fasta$|fastq$"  
then
grep ">" $1 | \
awk '{ 
	gsub(/>/,"",$1)
	split($1,a,"/")
	print $1 "\t" a[1] "\t" a[2] "\t" a[3]
}' >> $2
fi


