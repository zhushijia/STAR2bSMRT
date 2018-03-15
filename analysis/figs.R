library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1.txt',sep='\t')

folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/"
setwd(folder)
samples= dir()
samples = samples[samples!="combinedSamples"]

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


plotHeatmap = function(cc_fracs,ccExp)
{
	ccExp = log(ccExp+1)
	x = cbind( (1-cc_fracs) , 0 , 1-ccExp/max(ccExp) )
	x = x[order(x[,2]),]
	heatmap( x ,Rowv=NA,Colv=NA , scale="none")
}

setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/result/pipeline_nonAdjustNCjunc")
pdf("case.heatmap.pdf")
plotHeatmap( case_fracs , caseExp )
dev.off()

pdf("cont.heatmap.pdf")
plotHeatmap( cont_fracs , contExp )
dev.off()


plotExonCoexpress = function(cc_fracs, metrics=c("correlation","Euclidean"))
{
colnames(cc_fracs)=1:ncol(cc_fracs)
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
pdf("case.exonCoexpress.pdf")
plotExonCoexpress( case_fracs , "Euclidean" )
dev.off()

pdf("cont.exonCoexpress.pdf")
plotExonCoexpress( cont_fracs , "Euclidean" )
dev.off()


