echo "Input: $1"
echo "Output: $2"

echo -e "count\tchr\tstart\tend\tmotif" > $2

if echo $1 | grep "bam"; then
	samtools view $1 | \
	awk '{
	 n=split($19,a,",")
	 m=split($18,b,",")
	 if(a[2]!=-1)
	  for (i=1; i<=(n-1)/2; i++)
	   print $3 "\t" a[2*i] "\t"  a[2*i+1] "\t" b[i+1]
	}' | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' >> $2
fi

if echo $1 | grep "sam"; then 
	cat $1 | \
	awk '{
	 n=split($19,a,",")
	 m=split($18,b,",")
	 if(a[2]!=-1)
	  for (i=1; i<=(n-1)/2; i++)
	   print $3 "\t" a[2*i] "\t"  a[2*i+1] "\t" b[i+1]
	}' | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' >> $2

fi
