library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone2/setup")

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


hiPSC = unionGff(gffs[[1]],gffs[[2]])
humanBrain = unionGff(unionGff(gffs[[6]],gffs[[7]]),gffs[[8]])
fetal = gffs[[9]]
patient = unionGff(gffs[[3]],gffs[[4]])
mergeGffs = list( hiPSC=hiPSC , adult=humanBrain , fetal=fetal , patient=patient )





setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/result/pipeline_nonAdjustNCjunc")

ss = length(gffs)
png('gff_overlap.png',h=4200,w=4000)
par(mfrow=c(ss,ss))

cc = matrix(nrow=ss,ncol=ss,data=0)
for(i in 1:ss)
{
  for(j in 1:ss)
  {
    index = matchGff(gffs[[i]],gffs[[j]])
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



vennDiagramGff(gffs[c(1:5)] , filename=paste0("vennDiagram_",0,".png") )

for(i in 1:5)
{
	range = c(i,6:9)
	vennDiagramGff(gffs[range] , filename =paste0("vennDiagram_",samples[i],".png") )
	
}

vennDiagramGff(mergeGffs[1:3] , filename ="vennDiagram_hiPSC_humanbrain.png" )
vennDiagramGff(mergeGffs[2:4] , filename ="vennDiagram_patient_humanbrain.png" )









