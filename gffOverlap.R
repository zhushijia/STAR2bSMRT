samples= c("SZ_581","SZ_641","CT_2607","CT_553")
gff = list()
gff[[1]] = readGff( "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/581/Exp/isoform_ts2_td4.gff" , chrom='chr2' , s=50149082 , e=51255411 )
gff[[2]] = readGff( "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/641/Exp/isoform_ts4_td9.gff" , chrom='chr2' , s=50149082 , e=51255411 )
gff[[3]] = readGff( "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/2607/Exp/isoform_ts34_td7.gff" , chrom='chr2' , s=50149082 , e=51255411 )
gff[[4]] = readGff( "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/553/Exp/isoform_ts15_td23.gff" , chrom='chr2' , s=50149082 , e=51255411 )


exp = lapply( gff , function(x) {
  sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
}   )



setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/result")
png('gff_overlap.png',res=300,h=4200,w=4000)
par(mfrow=c(4,4))

cc = matrix(nrow=4,ncol=4,data=0)
for(i in 1:4)
{
  for(j in 1:4)
  {
    if(i!=j)
    {
      index = match(gff[[i]],gff[[j]])
      ei = exp[[i]][ !is.na(index) ] 
      ej = exp[[j]][ index[!is.na(index)] ]
      cc[i,j] = cor.test( ei , ej )$estimate
      plot(ei,ej,xlab=samples[i],ylab=samples[j],main=signif(cc[i,j],4) )
      
    } else {
      plot(ei,ej,col='white',xlab=samples[i],ylab=samples[j])
      
    } 
    
  }
}

dev.off()


library(VennDiagram)
tlist = lapply( gff , function(y) lapply( y , function(x) paste(paste(x[,1],x[,2])))   )
tall = unique( do.call(c, tlist ))
xx <- lapply( tlist , function(x) which(x%in%tall) )
names(xx) = samples
venn.diagram(xx, filename ="vennDiagram.png", height = 1800, width = 1800, fill=1:4 , alpha=0.3 ,imagetype = "png")


