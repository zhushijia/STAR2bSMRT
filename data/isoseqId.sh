echo "Input: $1"
echo "Output: $2"

echo -e "movieName\tZMWnumber\tsubread_start_end" > $2

grep ">" $1 | \
awk '{ 
	gsub(/>/,"",$1)
	split($1,a,"/")
	print a[1] "\t" a[2] "\t" a[3]
}' >> $2

