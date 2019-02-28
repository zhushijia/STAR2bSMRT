
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


getPSI = function(ef)
{
  # /sc/orga/projects/schzrnas/Erin/NRXN_projects/NRXNRNASeq/NRXN1_data/
  library(STAR2bSMRT,lib.loc="/sc/orga/projects//schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")
  annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
  parent="/sc/orga/projects/schzrnas/Erin/NRXN_projects/NRXNRNASeq/NRXN1_data"
  chrom = "chr2"
  s = 50147488
  e = 51259537
  SoutputDir = paste0(parent,'/',ef,'/analysis')
  SRalignment = paste0(SoutputDir,"/",gsub("Sample_","",ef),".merged.bam")
  SJout = paste0(gsub("Sample_","",ef),".final.sj.out.tab")
  #SRjunc = getJuncBySam( SRalignment, SoutputDir, chrom=chrom, s=s, e=e )
  SRjunc = getJuncBySJout( SJout=SJout,outputDir=SoutputDir,chrom=chrom,s=s,e=e )
  junc = SRjunc[['chr2']]
  
  jstop = function(exon) { subset(annotation , Exon==exon)$fullStop-1 }
  jstart = function(exon) { subset(annotation , Exon==exon)$fullStart+1 }
  PSI = function( exonLeft , exonMid , exonRight )
  {
    LeftIn = with( junc, count[ end == jstop(exonLeft) & start == jstart(exonMid) ] )
    RightIn = with( junc, count[ end == jstop(exonMid) & start == jstart(exonRight) ] )
    Ex = with( junc, count[ end == jstop(exonLeft) & start == jstart(exonRight) ] )
    
    if(length(LeftIn)==0) LeftIn=0
    if(length(RightIn)==0) RightIn=0
    if (exonLeft=='exon22')
      Ex = Ex + with( junc, count[ end == jstop('exon20') & start == jstart(exonRight) ] )
    if(length(Ex)==0) Ex=0
    
    #psi = (LeftIn+RightIn)/(LeftIn+RightIn+Ex)
    #if(is.na(psi)) psi=0
    psi = Ex/(LeftIn+RightIn+Ex)
    psi
  }
  
  psi = c()
  psi['3a'] =  PSI(exonLeft='exon2', exonMid='exon3a', exonRight='exon4')
  psi['3b'] =  PSI(exonLeft='exon2', exonMid='exon3b', exonRight='exon4')
  psi['4'] =   PSI(exonLeft='exon3a', exonMid='exon4', exonRight='exon5')
  psi['5'] =   PSI(exonLeft='exon4', exonMid='exon5', exonRight='exon6')
  
  psi['7a'] =  PSI(exonLeft='exon6', exonMid='exon7a', exonRight='exon8')
  psi['7b'] =  PSI(exonLeft='exon6', exonMid='exon7b', exonRight='exon8')
  
  psi['12'] =  PSI(exonLeft='exon11', exonMid='exon12', exonRight='exon13')
  psi['17'] =  PSI(exonLeft='exon16', exonMid='exon17', exonRight='exon18')
  psi['21'] =  PSI(exonLeft='exon20', exonMid='exon21', exonRight='exon22')
  psi['23a'] = PSI(exonLeft='exon22', exonMid='exon23a', exonRight='exon24')
  psi['23b'] = PSI(exonLeft='exon22', exonMid='exon23b', exonRight='exon24')
  
  psi
  
}


setwd("/sc/orga/projects/schzrnas/Erin/NRXN_projects/NRXNRNASeq/NRXN1_data")
samples = dir()[grep("Sample",dir())]
ratio = lapply(samples,function(sample) {
  cat(sample,'\n')
  getPSI(sample)
} )
ratio = do.call(rbind,ratio)

exp = lapply( samples, function(sample) {
  cat(sample,'\n')
 read.table( paste0(sample,"/analysis/",gsub("Sample_","",sample),".counts.txt") , header=T)
} )
geneExp = do.call(cbind, lapply(exp,function(x) x[,7]) )
ensembl_id = sapply( strsplit( as.character(exp[[1]][,1]),"[.]"), function(x) x[1] )
rownames(geneExp) = ensembl_id
colnames(geneExp) = samples

geneExp = data.frame( ensembl_id=ensembl_id, geneExp )
geneExp = merge(mapping,geneExp)
rownames(geneExp) = geneExp$symbol_id
geneExp = geneExp[,-c(1:4)]

condition = rep(1,length(samples))
condition[grep('NPC',samples)] = 2
outInfo = DESeq2.callDE(  geneExp , condition ) 

geneSd = apply(geneExp,2,sd)
geneExp = geneExp[,geneSd>0]

normExp = outInfo[,45:88]
nrxn = as.numeric( normExp[ grep("NRXN1",rownames(normExp)) , ] )
resExp = t( apply(normExp,1,function(x) residuals(lm(log(x+1)~log(nrxn+1))) ))

i = 8

p0 = apply(normExp,1,function(x) summary(lm( log(x+1)~ratio[,i] ))$coef[2,4]   )
adjp0 = p.adjust(p0)
paste( rownames(normExp)[order(p0)][1:1000] , collapse=" ")

p01 = apply(normExp,1,function(x) summary(lm( log(x+1)~nrxn ))$coef[2,4]   )
paste( rownames(normExp)[order(p01)][1:1000] , collapse=" ")

p1 = apply(normExp,1,function(x) summary(lm( log(x+1)~log(nrxn+1)+ratio[,i] ))$coef[3,4]   )
adjp1 = p.adjust(p1)
paste( rownames(normExp)[order(p1)][1:1000] , collapse=" ")


p2 = apply(resExp,1,function(x) summary(lm( x~ratio[,i] ))$coef[2,4]   )
adjp2 = p.adjust(p2)


ratio2 = ratio[,apply(ratio,2,sum,na.rm=T)>0]

pp1 = lapply( 1:ncol(ratio2) , function(i)
{
  cat(i,"\n")
  apply(normExp,1,function(x) summary(lm( log(x+1)~log(nrxn+1)+ratio2[,i] ))$coef[3,4]   )
} )


geneid = rownames(normExp)

enrich = lapply( 1:ncol(ratio2) , function(i)
{
  cat(i,"\n")
  gl = geneid[order(pp1[[i]])][1:500]
  enrichment = sapply( kegg , function(gs) { 
      p = fisher( gl, gs, geneid )$pvalue
      -log10(p)
    } )
  names(enrichment) = names(kegg)
  enrichment
} )

do.call(enrich)

xx= sapply(juncs,function(x)sum(x$count))
