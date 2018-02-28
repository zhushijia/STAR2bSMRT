library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

folder="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipeline_nonAdjustNCjunc/"
setwd(folder)
samples= dir()

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





setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/result/pipeline_nonAdjustNCjunc")

ss = length(gffs)
png('gff_overlap.png',h=4200,w=4000)
par(mfrow=c(ss,ss))

cc = matrix(nrow=ss,ncol=ss,data=0)
for(i in 1:ss)
{
  for(j in 1:ss)
  {
    index = match(gffs[[i]],gffs[[j]])
    ei = exp[[i]][ !is.na(index) ] 
    ej = exp[[j]][ index[!is.na(index)] ]
    cc[i,j] = cor.test( ei , ej )$estimate
    if(i!=j)
    {
      plot(ei,ej,xlab=samples[i],ylab=samples[j],main=signif(cc[i,j],4) )
    } else {
      plot(ei,ej,col='white',xlab=samples[i],ylab=samples[j])
      
    } 
    
  }
}

dev.off()



library(VennDiagram)
range = c(1:5)
tlist = lapply( gffs[range] , function(y) lapply( y , function(x) paste(paste(x[,1],x[,2])))   )
tall = unique( do.call(c, tlist ))
xx <- lapply( tlist , function(x) which(x%in%tall) )
names(xx) = samples[range]
venn.diagram(xx, filename =paste0("vennDiagram_",0,".png"), height = 1800, width = 1800, fill=1:5 , alpha=0.3 ,imagetype = "png")

for(i in 1:5)
{
	range = c(i,6:9)
	tlist = lapply( gffs[range] , function(y) lapply( y , function(x) paste(paste(x[,1],x[,2])))   )
	tall = unique( do.call(c, tlist ))
	xx <- lapply( tlist , function(x) which(x%in%tall) )
	names(xx) = samples[range]
	venn.diagram(xx, filename =paste0("vennDiagram_",samples[i],".png"), height = 1800, width = 1800, fill=1:5 , alpha=0.3 ,imagetype = "png")
}

