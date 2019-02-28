module load subread
gtf="/hpc/users/xzhus01/schzrnas/sjzhu/RNAseq/Reference/gencode/GRCh37.liftover.from.GRCh38.p10.Release26/annotation/gencode.v26lift37.annotation.gtf"
cd /sc/orga/projects/schzrnas/Erin/NRXN_projects/NRXNRNASeq/NRXN1_data
samples=$( ls | grep "Sample_" | sed s/Sample_//g )

for sample in $samples; do
bam="/sc/orga/projects/schzrnas/Erin/NRXN_projects/NRXNRNASeq/NRXN1_data/Sample_$sample/analysis/$sample.final.bam"
outputDir="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Erin_wholeGenome/Sample_$sample"
mkdir $outputDir
cd $outputDir
featureCounts -T 10 -t exon -f -O -p -a $gtf -o Sample_${sample}.exon.txt $bam
done





load("/hpc/users/xzhus01/schzrnas/sjzhu/Project/GWAS/GIGSEA/packageTest/GIGSEA/data/MSigDB.KEGG.Pathway.rda")
net = MSigDB.KEGG.Pathway$net
annot = MSigDB.KEGG.Pathway$annot
kegg = apply( net , 2, function(x) rownames(net)[x==1] )

fisher<-function( set1 , set2 , wholeset ) 
{
  k = sum(set1%in%set2)
  n = length(set1)    
  M = length(set2)  
  N = length(wholeset)
  d <- data.frame( in.interest = c(k, n-k) , not.interest = c(M-k, N-M-n+k) ) 
  row.names(d) <- c("In_category", "Not_in_category") 
  ft = fisher.test(d)     
  foldEnrichment = (k/n)/(M/N)         
  data.frame(k=k,n=n,M=M,N=N,pvalue=ft$p.value,foldEnrich=foldEnrichment)
}


library(org.Hs.eg.db)
symbol2entrez  = toTable(org.Hs.egSYMBOL)   # Map between Entrez Gene Identifiers and Gene Symbols
ensembl2entrez = toTable(org.Hs.egENSEMBL)  # Map Ensembl gene accession numbers with Entrez Gene identifiers
range = match( symbol2entrez$gene_id , ensembl2entrez$gene_id )
mapping = data.frame( symbol2entrez , ensembl2entrez[range,]  )
colnames(mapping) = c( 'entrez_id','symbol_id','entrez_id.1','ensembl_id')

require("DESeq2")

DESeq2.callDE <- function(  countData , condition , covariates=NULL , countThres=2 ) 
{
  # countData: count matrix with rowname being gene ids
  # condition: experimental conditions or groups
  # dds$condition <- factor(dds$condition, levels=c("untreated","treated"))
  # covariates should be a data.frame with column names representing variables 
  # countThres is threshold for the total read counts of one gene of all samples
  
  countData       <- countData[ rowSums(countData) >= countThres , ]  # filtering 
  
  if(is.null(covariates)) 
  {
    colData     <- data.frame( row.names=colnames(countData), condition )
    design      <- as.formula('~condition')
  } else {
    colData     <- data.frame( row.names=colnames(countData), condition , covariates )
    design      <- as.formula( paste0('~',paste(colnames(covariates),collapse='+'),'+condition') )
  }
  
  cat('The design is ' , as.character(design) , '\n')
  cat( nrow(countData) , 'genes survive the filtering\n')
  
  dds             <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=design )
  dds             <- DESeq(dds)
  res             <- results(dds)
  
  comparison      <- paste( levels(dds$condition)[2] , 'VS' , levels(dds$condition)[1] , sep='_')
  normalized.data <- counts(dds,normalized=TRUE)  # divided by sizeFactors(dds)
  log2fc          <- res$log2FoldChange
  lfcSE			<- res$lfcSE
  pvalue          <- res$pvalue
  fdr             <- res$padj
  
  outInfo <- data.frame( countData, normalized.data , log2fc , lfcSE , pvalue , fdr )
  colnames(outInfo) = c( paste( 'original',condition,sep='.' ) , paste( 'normalized',condition,sep='.') , paste(comparison,'log2fc',sep='_'), 'lfcSE' , 'pvalue' , 'fdr' ) 
  
  cat( sum(fdr<0.01,na.rm=T) , 'genes give fdr<0.01\n' )
  
  outInfo
}




annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)

library(foreach)
library(doMC)
registerDoMC(40)

setwd("/sc/orga/projects/schzrnas/Erin/NRXN_projects/NRXNRNASeq/NRXN1_data")
samples = dir()[grep("Sample",dir())]

exp1 = foreach( i=1:length(samples) ) %dopar% 
{
  sample = samples[i]
  cat(sample,'\n')
  read.table( paste0(sample,"/analysis/",gsub("Sample_","",sample),".counts.txt") , header=T)
}
geneExp = do.call(cbind, lapply(exp1,function(x) x[,7]) )
ensembl_id = sapply( strsplit( as.character(exp1[[1]][,1]),"[.]"), function(x) x[1] )
geneExp = data.frame( ensembl_id=ensembl_id, geneExp )
geneExp = merge(mapping,geneExp)
rownames(geneExp) = geneExp$symbol_id
geneExp = geneExp[,-c(1:4)]
colnames(geneExp) = samples

exp2 = foreach( i=1:length(samples) ) %dopar% 
{
  sample = samples[i]
  cat(sample,'\n')
  s = 50147488
  e = 51259537
  x = read.table( paste0("/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/data/Erin_wholeGenome/",
                         sample,"/",sample,".exon.txt") ,sep="\t", header=T)
  y = subset(x, Chr=='chr2' & Start<e & End>s ) 
  ind = lapply( 1:nrow(annotation) , function(i)
    which( y$Start==annotation$Stop[i] & y$End==annotation$Start[i] ) )
  ind2 = do.call(c, lapply(ind, function(a) a[1] ) )
  res = cbind(annotation,y[ind2,])[!is.na(ind2),]
  res
}
exonExp = do.call(cbind, lapply(exp2,function(x) x[,13]) )
colnames(exonExp) = samples
rownames(exonExp) = as.character(exp2[[1]][,1])

Exp = rbind(geneExp,exonExp)
condition = rep(1,length(samples))
condition[grep('NPC',samples)] = 2

outInfo = DESeq2.callDE(  Exp , condition ) 
colnames(outInfo)[1:44] = paste0('count.',samples)
colnames(outInfo)[45:88] = paste0('normalized.',samples)

normExp = outInfo[,45:88]
exon = normExp[ grepl('exon',rownames(normExp)) , ]
gene = normExp[ !grepl('exon',rownames(normExp)) , ]
nrxn = as.numeric( normExp[ grep("NRXN1",rownames(normExp)) , ] )
#exonRatio = t( apply(exon,1,function(x) x/nrxn ) )

setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Erin_wholeGenome/results")
write.table(outInfo,"normExp.txt",col.names=T,row.names=F,sep="\t",quote=F)

pp0 = apply(gene,1,function(x) summary(lm( x~nrxn ))$coef[2,4]   )
adjpp0 = p.adjust(pp0)
paste( rownames(gene)[adjpp0<1e-2] , collapse=" ")

pp1 = lapply( 1:nrow(exon) , function(i){
  cat(i,'\n')
  apply(gene,1,function(x) summary(lm( x~nrxn+as.numeric(exon[i,])+condition ))$coef[3,4]   )
})

adjpp1 = lapply(pp1,function(x) p.adjust(x) )
sg = lapply(adjpp1,function(fdr)  rownames(gene)[fdr<0.05] )
allSigGenes_fdr005 = unique( do.call(c,sg) ) # for GSEA


for (i in 1:nrow(exon))
{
  cat(i,'\n')
  info = data.frame( gene=rownames(gene), pval=pp1[[i]], fdr=adjpp1[[i]]  )
  info = info[order(info$pval),]
  setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/Erin_wholeGenome/results")
  fileName = paste0( 'parCorr_',rownames(exon)[i],'_withGeneExp.txt')
  write.table(info,fileName, sep='\t',col.names=T,row.names=F,quote=F)
}



i=20
p0 = apply(gene,1,function(x) summary(lm( x~nrxn+as.numeric(exon[i,]) ))$coef[3,4]   )
adjp1 = p.adjust(p1)
paste( rownames(gene)[p1<1e-2] , collapse=" ")


i=20
p1 = apply(gene,1,function(x) summary(lm( x~nrxn+as.numeric(exon[i,])+condition ))$coef[3,4]   )



i=20
p2 = apply(gene,1,function(x) summary(lm( x~as.numeric(exonRatio[i,])+condition ))$coef[2,4]   )
adjp2 = p.adjust(p2)
paste( rownames(gene)[p2<1e-2] , collapse=" ")


allgenes = rownames(gene)[ rownames(gene) %in% do.call(c,kegg) ]
gl = rownames(gene)[order(p1)][1:1000]
gl = gl[gl %in% allgenes]
enrichment = sapply( kegg , function(gs) { 
  p = fisher( gl, gs, allgenes )$pvalue
  p
} )
names(enrichment) = names(kegg)
#sort(enrichment,decreasing=T)
sort(p.adjust(enrichment),decreasing=T)




enrich = lapply( 1:length(pp1) , function(i)
{
  cat(i,"\n")
  allgenes = rownames(gene)
  #gl = rownames(gene)[order(pp1[[i]])][1:1000]
  gl = rownames(gene)[ p.adjust(pp1[[i]])<5e-2 ]
  enrichment = sapply( kegg , function(gs) { 
    p = fisher( gl, gs, allgenes )$pvalue
    p
  } )
  names(enrichment) = names(kegg)
  enrichment
} )



do.call(enrich)

xx= sapply(juncs,function(x)sum(x$count))
