library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

#annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1.txt',sep='\t')
annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1ExonAnnotations.txt',sep='\t',header=T)
annotation[2,2] = annotation[3,3]+1
annotation[25,3] = annotation[24,2]-1


folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/"
setwd(folder)
samples= dir()
samples= c("2607","553","581","641","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal")

gffs = lapply( samples , function(sample)
{
	files = dir( paste0(sample,"/Exp") )
	gff = paste0(sample,"/Exp/",files[grep("gff",files)])
	readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(gffs) = samples

exp = lapply( gffs , function(x) {
  sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
}   )

caseGff = unionGff(gffs[[3]],gffs[[4]])
contGff = unionGff(gffs[[1]],gffs[[2]])


case_sites = lapply(caseGff,function(y) apply(y,1,function(z) z[2]:z[3] ) )
cont_sites = lapply(contGff,function(y) apply(y,1,function(z) z[2]:z[3] ) )
annot_sites = apply(annotation,1,function(z) z[2]:z[3] )

caseExp = sapply( 3:4 , function(i) {
	ee = exp[[i]][ matchGff( caseGff , gffs[[i]] ) ]
	ee[is.na(ee)] = 0
	ee
} )
contExp = sapply( 1:2 , function(i) {
	ee = exp[[i]][ matchGff( contGff , gffs[[i]] ) ]
	ee[is.na(ee)] = 0
	ee
} )


frac = function(gs,annot_sites)
{
	sapply(annot_sites,function(x) {
		max( sapply(gs,function(y) mean(x%in%y)) )
	} )
}

case_fracs = do.call(rbind,lapply( case_sites , function(gs) frac(gs,annot_sites) ) )
cont_fracs = do.call(rbind,lapply( cont_sites , function(gs) frac(gs,annot_sites) ) )
colnames(case_fracs) = as.character(annotation$Exon)
colnames(cont_fracs) = as.character(annotation$Exon)



if(0)
{
	library(Biostrings)
	ref = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
	genome = readDNAStringSet(ref)
	caseSeq = generateSeq( genome , isoform=caseGff )
	contSeq = generateSeq( genome , isoform=contGff )

	as.numeric(case_fracs[caseSeq$translated,25])
	as.numeric(case_fracs[caseSeq$translated,2])
}


plotHeatmap = function(cc_fracs,ccExp)
{
	ccExp = log(ccExp+1)
	x = cbind( (1-cc_fracs) , 0 , 1-ccExp/max(ccExp) )
	x = x[order(x[,2]),]
	heatmap( x ,Rowv=NA,Colv=NA , scale="none")
}

setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/result/pipeline_nonAdjustNCjunc")
pdf("case.heatmap2.pdf")
plotHeatmap( case_fracs , caseExp )
dev.off()

pdf("cont.heatmap2.pdf")
plotHeatmap( cont_fracs , contExp )
dev.off()


plotExonCoexpress = function(cc_fracs, metrics=c("correlation","Euclidean"))
{
#colnames(cc_fracs)=1:ncol(cc_fracs)
if( metrics=="correlation")
{
	sds = apply(cc_fracs,2,sd)
	cc_fracs = cc_fracs[,sds>0]
	heatmap( 1-cor(cc_fracs) ,Rowv=NA,Colv=NA , scale="none")
}
if( metrics=="Euclidean")
	heatmap( as.matrix(dist(t(cc_fracs))) ,Rowv=NA,Colv=NA , scale="none")
}

setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/result/pipeline_nonAdjustNCjunc")
pdf("case.exonCoexpress2.pdf")
plotExonCoexpress( case_fracs , "Euclidean" )
dev.off()

pdf("cont.exonCoexpress2.pdf")
plotExonCoexpress( cont_fracs , "Euclidean" )
dev.off()


