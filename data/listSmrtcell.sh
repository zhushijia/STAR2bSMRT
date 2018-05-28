echo "Option: $1"
echo "Input: $2"
echo "Output: $3"

if [ $1 == "-uniq" ]; then
echo "... output unique smrtcells and counts"
echo -e "count\tsmrtcell" > $3
samtools view $2 | \
awk '{
n=split($1,a,"/")
print a[1]
}' | sort | uniq -c >> $3
fi


if [ $1 == "-all" ]; then
echo "... output all smrtcells and read length"
echo -e "smrtcell\treadlen" > $3
samtools view $2 | \
awk '{
n=split($1,a,"/")
print a[1] "\t" length($10)
}' >> $3
fi

