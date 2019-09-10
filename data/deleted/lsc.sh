lsc() {

local LR=$1
local SR1=$2
local SR2=$3
local outputDir=$4

cd $outputDir
module load bowtie2
/hpc/users/zhus02/setup/LSC-2.0/bin/runLSC.py --long_reads $LR --short_reads $SR1 $SR2 --threads 10 --short_read_file_type fq --output lsc

}

