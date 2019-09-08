
kallistoMap ()
{
	local idx=$1
	local R1=$2
	local R2=$3
	local outputDir=$4

	mkdir -p $outputDir
	cd $outputDir
	module load kallisto
	#kallisto index -i transcripts.idx $1
	kallisto quant -i $idx -o output -b 10 $R1 $R2
}


########################################################################################################################
##########################################   Hiseq  ######################################################################
########################################################################################################################

parent="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/NRXN1_abgRef/kallisto"
idx="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/NRXN1_abgRef/transcripts.idx"

SR1="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_581-2-1-FBN/Sample_581-2-1-FBN.R1.fastq"
SR2="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_581-2-1-FBN/Sample_581-2-1-FBN.R2.fastq"
outputDir=$parent/581-2-1
kallistoMap $idx $SR1 $SR2 $outputDir


SR1="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_641-6-2-FBN/Sample_641-6-2-FBN.R1.fastq"
SR2="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_641-6-2-FBN/Sample_641-6-2-FBN.R2.fastq"
outputDir=$parent/641-6-2
kallistoMap $idx $SR1 $SR2 $outputDir


SR1="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_2607-2-1-FBN/Sample_2607-2-1-FBN.R1.fastq"
SR2="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_2607-2-1-FBN/Sample_2607-2-1-FBN.R2.fastq"
outputDir=$parent/2607-2-1
kallistoMap $idx $SR1 $SR2 $outputDir


LR="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24459_553/polished_high_qv_consensus_isoforms.fasta"
SR1="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_553-1-1-FBN/Sample_553-1-1-FBN.R1.fastq"
SR2="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/ShortReads/Sample_553-1-1-FBN/Sample_553-1-1-FBN.R1.fastq"
outputDir=$parent/553-1-1
kallistoMap $idx $SR1 $SR2 $outputDir



