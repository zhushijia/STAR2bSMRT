library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")
library(Biostrings)

gff2tag = function(gff)
{
  lapply( gff,function(x) paste( x$chr[1] , paste( paste( x$end[-nrow(x)] , x$start[-1]  ),collapse="; "),sep=": " ) )
}

get_tofuInfo = function(LRphqv, TOFUgroup, TOFUgff)
{
	
  hq = readDNAStringSet(LRphqv)
	seq_name = sapply( strsplit(names(hq)," "), function(x) x[1])
	#seq_name = paste0( jobid , '_' , seq_name)
	seq_info = sapply( strsplit(names(hq)," "), function(x) x[2])
	seq_info_list = strsplit(seq_info,";")
	isoform                  = sapply( seq_info_list , function(x) gsub('isoform=','',x[1]) )
	full_length_coverage     = sapply( seq_info_list , function(x) as.integer(gsub('full_length_coverage=','',x[2])) )
	non_full_length_coverage = sapply( seq_info_list , function(x) as.integer(gsub('non_full_length_coverage=','',x[3])) )
	isoform_length           = sapply( seq_info_list , function(x) as.integer(gsub('isoform_length=','',x[4])) )
	hq_info = data.frame(seq_name , isoform , full_length_coverage , non_full_length_coverage , isoform_length)
	
	tofu = read.table(TOFUgroup,header=F,sep='\t')
	tofu_group = strsplit( as.character(tofu[,2]) , ',' )
	names(tofu_group) = as.character(tofu[,1])
	exp = do.call( c , lapply( tofu_group, function(x) 
	  sum(subset(hq_info,seq_name%in%x)$full_length_coverage) ) )
	
	tofu_info = data.frame(transcript=names(exp),full_length_coverage=exp)
	tofu_gff = readGff( TOFUgff , chrom='chr2' , s=50149082 , e=51255411 )
	names(tofu_gff) = sapply( strsplit(names(tofu_gff)," "), function(x) gsub(";","",x[4]))
	tofu_gff = tofu_gff[ match( as.character(tofu_info$transcript), names(tofu_gff) ) ]
	tofu_exp = tofu_info$full_length_coverage
	all(names(tofu_gff)==as.character(tofu_info[,1]))
	tofu_junc = gff2tag( tofu_gff )
	
	list(gff=tofu_gff,exp=tofu_exp,junc=tofu_junc)
	
}

gabaLRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/028024/polished_high_qv_consensus_isoforms.fasta"
gabaTOFUgroup = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/641_GABA/LR/tofu/tofu.collapsed.group.txt"
gabaTOFUgff = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/641_GABA/LR/tofu/tofu.collapsed.gff"

ngn2LRphqv="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Run6/027959/polished_high_qv_consensus_isoforms.fasta"
ngn2TOFUgroup = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/641_NGN2/LR/tofu/tofu.collapsed.group.txt"
ngn2TOFUgff = "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/641_NGN2/LR/tofu/tofu.collapsed.gff"

gaba_tofuInfo = get_tofuInfo(gabaLRphqv, gabaTOFUgroup, gabaTOFUgff)
ngn2_tofuInfo = get_tofuInfo(ngn2LRphqv, ngn2TOFUgroup, ngn2TOFUgff)
    

annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
fbnGff = readGff( "/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/641/adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100/isoform_ts94_td0.gff" , 
                  chrom='chr2' , s=50149082 , e=51255411 )
translated = sapply( strsplit(names(fbnGff),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
fbnGff =  fbnGff[ translated ]
fbnExp = sapply( strsplit(names(fbnGff),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
fbnGff = fbnGff[fbnExp>=7]
fbnExp = fbnExp[fbnExp>=7]
fbnJunc = gff2tag(fbnGff)

mutant = sapply( fbnJunc , function(x) grepl("50149389 50318461",x) )
gabaExp = gaba_tofuInfo$exp[ match(fbnJunc, gaba_tofuInfo$junc) ] 
ngn2Exp = ngn2_tofuInfo$exp[ match(fbnJunc, ngn2_tofuInfo$junc) ]
xx = data.frame( fbn=fbnExp , gaba=gabaExp , ngn2=ngn2Exp )
yy = log10(xx+1)          


setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/comparisonResult/adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100")
pdf('FBN_GAGA_NGN2_641.pdf')
par(mfrow=c(2,2))
plot( yy$fbn , yy$gaba , pch=16, xlab="log10 FBN read counts", ylab="log10 GABA read counts")
abline(lm(yy$gaba~yy$fbn))
points( yy$fbn[mutant] , yy$gaba[mutant] , pch=16, col='red' )

plot( yy$fbn , yy$ngn2 , pch=16, xlab="log10 FBN read counts", ylab="log10 NGN2 read counts")
abline( lm(yy$ngn2~yy$fbn) )
points( yy$fbn[mutant] , yy$ngn2[mutant] , pch=16, col='red' )
dev.off()

nrow(yy)
sum(mutant)

sum( !is.na(yy$gaba) )
sum( !is.na(yy$gaba) & mutant )

sum( !is.na(yy$ngn2) )
sum( !is.na(yy$ngn2) & mutant )

cor.test( yy$fbn , yy$gaba )
cor.test( yy$fbn , yy$ngn2 )
