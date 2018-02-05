echo "Input: $1"
echo "Output: $2"

if echo $1 | grep "sam"; then
	grep "jI:B:i" $1 | cut -f1,3,4,19 | sed s/jI:B:i,//g  > $2.txt
	samtools view -bS $1 | bamToBed -i > $2.bed
	echo -e "id\tchr\tstrand\tstart\tend\tjunc" > $2
	paste $2.bed $2.txt | awk '{print $4 "\t" $1 "\t" $6 "\t" $2 "\t" $3"\t" $10 }' >> $2
	rm $2.txt
	rm $2.bed
fi

if echo $1 | grep "bam"; then
	samtools view $1 | grep "jI:B:i" $1 | cut -f1,3,4,19 | sed s/jI:B:i,//g  > $2.txt
	bamToBed -i $1 > $2.bed
	echo -e "id\tchr\tstrand\tstart\tend\tjunc" > $2
	paste $2.bed $2.txt | awk '{print $4 "\t" $1 "\t" $6 "\t" $2 "\t" $3"\t" $10 }' >> $2
	rm $2.txt
	rm $2.bed
fi



