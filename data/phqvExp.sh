echo "Input: $1"
echo "Output: $2"

echo -e "id\tfull_length_coverage\tnon_full_length_coverage\tisoform_length" > $2

grep ">" $1 | \
awk '{ 
	split($2,a,";")
	#gsub(/isoform=/,"",a[1])
	gsub(/>/,"",$1)
	gsub(/full_length_coverage=/,"",a[2])
	gsub(/non_full_length_coverage=/,"",a[3])
	gsub(/isoform_length=/,"",a[4])
	print $1 "\t" a[2] "\t" a[3] "\t" a[4]
}' >> $2

