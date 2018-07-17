library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")
annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)



parent="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping"
chrom = "chr2"
s = 50147488
e = 51259537
SoutputDir = paste0(parent,'/',ef)
SRalignment = paste0(SoutputDir,"/alignments.bam")
SRjunc = getJuncBySam( SRalignment, SoutputDir, chrom=chrom, s=s, e=e )
junc = SRjunc[['chr2']]

junc$count[ junc$end == (subset(annotation , Exon=='exon2')$fullStop-1) ]



