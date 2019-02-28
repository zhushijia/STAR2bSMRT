library(STAR2bSMRT,lib.loc="/hpc/users/xzhus01/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone4/setup")

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters"
parai1 = "adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_TRUE_fuzzyMatch_100"
parai2 = "adjustNCjunc_TRUE_fixedMatchedLS_FALSE_useSJout_FALSE_fuzzyMatch_100"
#samples = c( paste0( c("581_GABA_fbnSR","581_NGN2_fbnSR","641_NGN2"),"/",parai2 ), paste0(c("581","641"),"/",parai2) )

parai = parai2
samples = paste0( c("553_GABA","553_NGN2","553"),"/", parai )
sampleName = c("553_GABA","553_NGN2","553_FBN")


gffs = lapply( samples , function(sample)
{
    files = dir( paste0(folder,"/",sample,"/") )
    gff = paste0(folder,"/",sample,"/",files[grepl("gff$",files)])
    readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
} )
names(gffs) = sampleName

exp = lapply( gffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
}   )

translated = lapply( gffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
}   )

for( thres in 1:10 )
{
  translatedGffs = lapply( 1:length(gffs) , function(i) gffs[[i]][ translated[[i]] & exp[[i]]>=thres ]  )
  translatedExps = lapply( 1:length(gffs) , function(i) exp[[i]][ translated[[i]]  & exp[[i]]>=thres ] )
  names(translatedGffs) = sampleName
  names(translatedExps) = sampleName
  num = sapply(translatedGffs,length)
  
  setwd( paste0(folder,"/comparisonResult/",parai ) )
  vennDiagramGff(translatedGffs , filename=paste0("vennDiagram_553_GABA_GLUT_covThres",thres,".png") )
}





ind = function(x1,x2)
{
    tmp = matchGff(x1,x2)
    tmp[!is.na(tmp)]
}

indlen = function(i,j,n)
{
    length( ind( translatedGffs[[i]][ order(translatedExps[[i]],decreasing=T)[1:n] ] , 
                 translatedGffs[[j]][ order(translatedExps[[j]],decreasing=T)[1:n] ] ) )
}

indlen2 = function(j,n)
{
    length( ind( translatedGffs[[3]] , 
                 translatedGffs[[j]][ order(translatedExps[[j]],decreasing=T)[1:n] ] ) )
}

x1 = sapply( 1:min(num) , function(n) indlen(1,3,n)/n )
x2 = sapply( 1:min(num) , function(n) indlen(2,3,n)/n )
x3 = sapply( 1:min(num) , function(n) indlen(2,1,n)/n )
y1 = sapply( 1:min(num) , function(n) indlen2(1,n)/n )
y2 = sapply( 1:min(num) , function(n) indlen2(2,n)/n )


library(ggplot2)
DF <- data.frame(
    x = c( x1, x2, x3 ),
    y = rep( as.character( ceiling(c(1:min(num)/10) )) , 3),
    z = rep( 1:3, each=min(num)),
    stringsAsFactors = FALSE
)
setwd( paste0(folder,"/comparisonResult/",parai ) )
pdf("overlapBarplot_553_GABA_GLUT.pdf")
ggplot(DF, aes(y, x, fill=factor(z))) + geom_boxplot()

plot(NA,xlim=c(1,81),ylim=c(0,1) ) 
points( x1 , col='red' , pch=15 )
points( x2 , col='green' , pch=16 )
points( x3 , col='blue'  , pch=17)
lines( x1 , col='red' , lwd=2 )
lines( x2 , col='green' , lwd=2 )
lines( x3 , col='blue'  , lwd=2 )

dev.off()



setwd( paste0(folder,"/comparisonResult/",parai ) )

wilcox.test(x1,x2,paired=T,'greater')
wilcox.test(y1,y2,paired=T,'greater')
