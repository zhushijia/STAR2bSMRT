library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone3/setup")
library(Biostrings)
library(foreach)
library(doMC)

cores = 30
thresSR=c(1:100) 
thresDis=c(1:20)
adjustNCjunc=TRUE  # adjustNCjunc=TRUE/FALSE and fixedMatchedLS=FALSE ts=70 td=13 num=282 frac=0.497, when using SJ.out.tab
fixedMatchedLS=FALSE
chrom="chr17"
s= 90036900
e = 91089605
registerDoMC(cores)


outputDir="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/STAR2bSMRT"
LoutputDir = paste0(outputDir,"/LR")
SoutputDir = paste0(outputDir,"/SR")
EoutputDir = paste0(outputDir,"/Exp")

system( paste0( "mkdir -p " , LoutputDir ) )
system( paste0( "mkdir -p " , SoutputDir ) )
system( paste0( "mkdir -p " , EoutputDir ) )

SRalignment = paste0(SoutputDir,"/alignments.bam")
LRalignment = paste0(LoutputDir,"/Aligned.out.sam")


SRjunc = getJuncBySam( SRalignment , SoutputDir , chrom="chr17" , s= 90036900 , e = 91089605 )
SRjunc = getJuncBySJout( SJout="SJ.out.tab" , SoutputDir , chrom="chr17" , s= 90036900 , e = 91089605 )

SJ.out.tab = read.table( paste0(SoutputDir,"/SJ.out.tab") , sep="\t")
SJ.out.tab = with( SJ.out.tab , data.frame(count=V7, chr=V1, start=V2, end=V3, motif=V5) )
SRjunc$chr17 = subset(SJ.out.tab , chr=="chr17" & start >= 90036900 & end <= 91089605 )


LRinfo = getLRinfo( LRalignment , NULL , LoutputDir , chrom="chr17" , s= 90036900 , e = 91089605 , jI=TRUE)
LRread = LRinfo$LRread
LRjunc = LRinfo$LRjunc
LRtag = LRinfo$LRtag

if( fixedMatchedLS )
{
  matchedLS = matchLSjunc( LRjunc , SRjunc )
} else {
  matchedLS = NULL
}

score = gridSearch( LRjunc , SRjunc , thresSR , thresDis , adjustNCjunc , matchedLS )

ij = which( score==max(score) , arr.ind=T )
ts = thresSR[ ij[1,1] ]
td = thresDis[ ij[1,2] ]
cat( ts , td , score[ij] , '\n ')

correction = generateCorrectedIsoform( LRjunc , SRjunc, LRtag , LRread  , ts , td , matchedLS )
correction[[1]]$frac
correction[[1]]$num

#######################################################################################################
setwd( EoutputDir )
pdf( paste0( "JuncExp_LR_ts",ts,"_td",td,".pdf") )

juncExp = do.call( rbind, lapply( correction , function(x) x$LSjuncCount ))
lrCount = log10(juncExp$lrCount)
srCount = log10(juncExp$srCount)
juncCorr = cor.test(srCount,lrCount,method='spearman')$estimate
cols = sapply( juncExp$motif , function(x) ifelse(x==0,1,2) )
cols[juncExp$motif==1] = 3
plot( lrCount , srCount , col=cols , pch=17 , main=paste0("JuncExp by Long and Short Reads: r=",signif(juncCorr,3)) ,  xlab="Log10 Long Read" , ylab="Log10 Short Read"  )
abline(lm( srCount~lrCount ))

par(mfrow=c(2,1))
log10fc = lrCount - srCount
JuncNames = paste(juncExp$start , juncExp$end)
barplot( log10fc , cex.names=0.6 , col=cols , ylab="log10(lrCount/srCount)", names=JuncNames , las=3 )

dev.off()

#######################################################################################################

gffName = paste0( "isoform_ts",ts,"_td",td,".gff")
exonList = juncToExon( juncList=correction[[chrom]]$isoform , s=90036900 , e=91089605 , exp=correction[[chrom]]$exp )
writeGff( isoform=exonList , file = gffName )


#######################################################################################################
# check consistency with human adult

#annotation = read.table("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/Mouse_NRXN1_Annotation.txt",sep="\t",header=T)
annotation = read.table("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Mouse_NRXN/Mouse_NRXN1_Annotation_matchedHuman.txt",sep="\t",header=T)
#annotation[6,2] = 90704367
mouseAnnot = annotation

annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1ExonAnnotations.txt',sep='\t',header=T)
annotation[2,2] = annotation[3,3]+1
annotation[25,3] = annotation[24,2]-1
humanAnnot = annotation

# load gff from the human analysis
humanAdultGff = unionGff(unionGff(gffs[[6]],gffs[[7]]),gffs[[8]])

#######################################################################################################
# check consistency with human adult

frac = function(gs,annot_sites)
{
  sapply(annot_sites,function(x) {
    max( sapply(gs,function(y) mean(x%in%y)) )
  } )
}

humanAnnot_sites = apply(humanAnnot,1,function(z) z[2]:z[3] )
human_sites = lapply(humanAdultGff,function(y) apply(y,1,function(z) z[2]:z[3] ) )
human_frac = do.call(rbind,lapply( human_sites , function(gs) frac(gs,humanAnnot_sites) ) )
colnames(human_frac) = as.character(humanAnnot$Exon)
#colnames(human_frac) = gsub("a|b$","",colnames(human_frac))
human_frac = human_frac[,1:25]

mouseGff = exonList
mouseAnnot_sites = apply(mouseAnnot,1,function(z) z[2]:z[3] )
mouse_sites = lapply(mouseGff,function(y) apply(y,1,function(z) z[2]:z[3] ) )
mouse_frac = do.call(rbind,lapply( mouse_sites , function(gs) frac(gs,mouseAnnot_sites) ) )
colnames(mouse_frac) = as.character(mouseAnnot$Exon)
#colnames(mouse_frac) = gsub("a|b$","",colnames(mouse_frac))


human_frac1 = round(human_frac)
mouse_frac1 = round(mouse_frac)

human_tag1 = apply(human_frac1,1,paste,collapse=",")
mouse_tag1 = apply(mouse_frac1,1,paste,collapse=",")

human_frac1 = human_frac1[order(human_tag1),]
mouse_frac1 = mouse_frac1[order(mouse_tag1),]

human_frac2 = t(apply(human_frac,1,function(x) tapply(x,gsub("a|b$","",colnames(human_frac)),function(y) round(max(y)) )))
mouse_frac2 = t(apply(mouse_frac,1,function(x) tapply(x,gsub("a|b$","",colnames(mouse_frac)),function(y) round(max(y)) )))
colnames(mouse_frac2) == colnames(mouse_frac2)

human_tag2 = apply(human_frac2,1,paste,collapse=",")
mouse_tag2 = apply(mouse_frac2,1,paste,collapse=",")

length(unique(human_tag2))
length(unique(mouse_tag2))
sum( unique(human_tag2) %in% unique(mouse_tag2) )



setwd( EoutputDir )
require('gplots')
pdf("hiPSC_mouse_isoform.pdf")
heatmap.2(-human_frac1 , main="human adult isoforms", col=heat.colors(256), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
heatmap.2(-mouse_frac1 , main="mouse isoforms", col=heat.colors(256), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none",keysize = 1.2,cexRow=0.5,dendrogram="none",Rowv=NA,Colv=NA)
dev.off()


#######################################################################################################
# following are potentially useful

annot_juncs = data.frame( js=annotation[2:nrow(annotation),2]+1 , je=annotation[1:(nrow(annotation)-1),1]-1 )
tag = lapply( exonList , function(x) apply(x,1,function(y) paste(y[2:3],collapse=",")) )
annot_tag = apply(annotation,1,function(y) paste(y[1:2],collapse=",") )

exonNames = lapply(tag,function(x) as.character(annotation$Annotation)[ match(x,annot_tag) ]   )
sapply( exonNames, function(x) all(!is.na(x)) )
table(do.call(c,sapply( 1:length(tag), function(i) tag[[i]][is.na(exonNames[[i]])] )))



hiPSCannot_tag = apply(humanAnnot,1,function(y) paste(y[2:3],collapse=",") )
hiPSCexonNames = lapply(hiPSCtag,function(x) as.character(humanAnnot$Exon)[ match(x,hiPSCannot_tag) ]   )
humanRef = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
genome = readDNAStringSet(humanRef)
humanSeq = generateSeq( genome , isoform=hiPSCGff )

