folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/"
samples= c("2607","553","581","641","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal")

parameters = dir( paste0(folder,"/fetal") )

for(parai in parameters)  
{
  
  cat(parai,"\n")
  
  outputDir = paste0(folder,"/comparisonResult","/",parai)
  system( paste( "mkdir -p" , outputDir ) )
  x = lapply( samples , function(sample) read.table( paste0(folder,"/",sample,"/",parai,"/summary.txt"),header=T) ) 
  x = do.call(rbind,x)
  y = data.frame( x , translatedFrac = x$translated/x$isoformNum )
  y = rbind( y , apply(y,2,mean) )
  info = data.frame( samples=c(samples,"mean") , signif(y,3) )
  
  outputDir = paste0(folder,"/comparisonResult","/",parai)
  system( paste( "mkdir -p" , outputDir ) )
  setwd(outputDir)
  write.table( info , "summary.csv" , sep="," , col.names=T,row.names=F,quote=F  )
  
}
