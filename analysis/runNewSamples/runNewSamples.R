library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
SoutputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping/EF50"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/EF50_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/EF50_R2_001.fastq.gz"
genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
starShort( genomeDir , SR1 , SR2 , SoutputDir )



library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")
SoutputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping/EF12"
SR1="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/EF12_R1_001.fastq.gz"
SR2="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/EF12_R2_001.fastq.gz"
genomeDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/IDPtest_ErinData/starShort/genomeDir_1pass"
starShort( genomeDir , SR1 , SR2 , SoutputDir )

