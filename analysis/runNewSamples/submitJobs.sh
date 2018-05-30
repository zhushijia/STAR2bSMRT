
parent="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035"
cd $parent
samples=$( ls | grep gz )
samples=$(for sample in $samples; do echo ${sample%_R*_*} ; done | sort | uniq )

jobName="starShort1"
mkdir $parent/$jobName

for sample in $samples ; do

genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
gtf="/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/gencode/GRCh37.liftover.from.GRCh38.p10.Release26/annotation/gencode.v26lift37.annotation.gtf"
R1="${parent}/${sample}_R1_001.fastq.gz"
R2="${parent}/${sample}_R2_001.fastq.gz"
outputDir="${parent}/mapping/${sample}"
mkdir -p $outputDir

echo "
#!/bin/bash
#BSUB -J $sample
#BSUB -W 8:00
#BSUB -n 10
#BSUB -q premium
#BSUB -P acc_schzrnas
#BSUB -o $outputDir/$jobName.o
#BSUB -e $outputDir/$jobName.e
#BSUB -u shijia.zhu@mssm.edu
#BSUB -R rusage[mem=5000]
#BSUB -R span[hosts=1]
#BSUB -L /bin/bash


. /hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/mappingCodes/star/starShort1.sh $genomeDir $R1 $R2 $outputDir
cd $outputDir
featureCounts -T 10 -t exon -f -O -p -a $gtf -o exon.txt alignments.bam
grep ENSG00000179915 exon.txt > NRXN1.exon.txt

" > $parent/$jobName/$sample.lsf

bsub < $parent/$jobName/$sample.lsf

done 


