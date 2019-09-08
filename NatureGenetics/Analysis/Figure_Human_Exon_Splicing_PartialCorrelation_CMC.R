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





require("DESeq2")

DESeq2.normalize <- function(  countData , countThres=2 ) 
{
    # countData: count matrix with rowname being gene ids
    # countThres is threshold for the total read counts of one gene of all samples
    
    countData       <- countData[ rowSums(countData) >= countThres , ]  # filtering 
    colData         <- data.frame(rep("null",ncol(countData)))
    dds             <- DESeqDataSetFromMatrix(countData=countData, colData=colData , design=as.formula('~1') )
    dds             <- DESeq(dds)
    normalized.data <- counts(dds,normalized=TRUE)  # divided by sizeFactors(dds)
    normalized.data
    
}


library(org.Hs.eg.db)
symbol2entrez  = toTable(org.Hs.egSYMBOL)   # Map between Entrez Gene Identifiers and Gene Symbols
ensembl2entrez = toTable(org.Hs.egENSEMBL)  # Map Ensembl gene accession numbers with Entrez Gene identifiers
range = match( symbol2entrez$gene_id , ensembl2entrez$gene_id )
mapping = data.frame( symbol2entrez , ensembl2entrez[range,]  )
colnames(mapping) = c( 'entrez_id','symbol_id','entrez_id.1','ensembl_id')



annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
nrxn1 = read.table('/hpc/users/xzhus01/schzrnas/sjzhu/RNAseq/Reference/GRCh38.84/NRXN1_hg38.txt',sep='\t')
exon = subset( nrxn1 , V3=='exon' )
exon_id = sapply( strsplit( as.character(exon[,9]), ';' ) , function(x) x[16] )
exon_id = gsub(' exon_id ' , '' , exon_id)
exon    = exon[grep('ENSE',exon_id), ]
exon_id = unique(exon_id[grep('ENSE',exon_id)])

library(foreach)
library(doMC)
registerDoMC(40)

p1 = "/sc/orga/projects/CommonMind/data/CMC_DATA_REDO_MAY_2018/HBCC_DLPFC"
p2 = "/sc/orga/projects/CommonMind/data/CMC_DATA_REDO_MAY_2018/MSSM_PENN_PITT_DLPFC"
f1 = dir(p1)[grep("CMC_HBCC",dir(p1))]
f2 = dir(p2)[grep("_RNA_",dir(p2))]
pf1 = paste0(p1,"/",f1)
pf2 = paste0(p2,"/",f2)
pf = c(pf1,pf2)

exp1 = foreach( i=1:length(pf) ) %dopar% 
{
  path = pf[i]
  sample = basename(pf[i])
  cat(sample,'\n')
  x = read.table( paste0(path,"/Processed/RAPiD.3_0_0/featureCounts/",
                     sample,".primary.txt") , header=T)
  x[,c(1,7)]
}

geneExp = do.call(cbind, lapply(exp1,function(x) x[,2]) )
colnames(geneExp) = basename(pf)
rownames(geneExp) = as.character(exp1[[1]][,1])
ensembl_id = sapply( strsplit( as.character(exp1[[1]][,1]),"[.]"), function(x) x[1] )
geneExp = data.frame( ensembl_id=ensembl_id, geneExp )
geneExp = merge(mapping,geneExp)
num = table(geneExp$symbol_id)
geneExp = geneExp[! geneExp$symbol_id %in% names(num[num>1]) , ]
rownames(geneExp) = geneExp$symbol_id
geneExp = geneExp[,-c(1:4)]
#colnames(geneExp) = basename(pf)

exp2 = foreach( i=1:length(pf) ) %dopar% 
{
    path = pf[i]
    sample = basename(pf[i])
    cat(sample,'\n')
    x = read.table( paste0(path,"/Processed/RAPiD.3_0_0/featureCounts/",
                           sample,".exon.txt") , header=T)
    range = grep("ENSG00000179915",as.character(x$Geneid))
    x[range,7]
}

exonExp = do.call(cbind, exp2 )
colnames(exonExp) = basename(pf)
rownames(exonExp) = paste0('exon',1:nrow(exonExp))

Exp = rbind(geneExp,exonExp)
setwd("/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/result/ExonExlucsionPartialCorr")
#save(Exp, file='CMC_Exp.rd')

normExp = DESeq2.normalize(  Exp , countThres=100 ) 



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
