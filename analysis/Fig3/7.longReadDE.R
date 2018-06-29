library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone4/setup")

folder="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters2/"
samples= c("2607","553","581","641","642","NRXN_adult_dlPFC1_12","adult_dlPFC1_10","adult_dlPFC1_13","fetal")
sampleNames = c("Control 1","Control 2","3' Del 1","3' Del 2","Control 3","Adult 1","Adult 2","Adult 3","Fetal")


parameters = dir( paste0(folder,"/fetal") )

for(parai in parameters)  
{
  
  cat(parai,"\n")
  
  gffs = lapply( samples , function(sample)
  {
    files = dir( paste0(folder,"/",sample,"/",parai) )
    gff = paste0(folder,"/",sample,"/",parai,"/",files[grepl("gff$",files)])
    readGff( gff , chrom='chr2' , s=50149082 , e=51255411 )
  } )
  names(gffs) = sampleNames
  
  translated = lapply( gffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
  }   )
  
  
  translatedGffs = lapply( 1:4 , function(i) gffs[[i]][ translated[[i]] ]  )
  translatedExps = lapply( translatedGffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
  }   )
  names(translatedGffs) = sampleNames[1:4]
  names(translatedExps) = sampleNames[1:4]
  
  contGff =  unionGff(translatedGffs[[1]],translatedGffs[[2]])
  caseGff =  unionGff(translatedGffs[[3]],translatedGffs[[4]])
  allGff =  unionGff(caseGff,contGff)
  
  expMatrix = sapply( 1:4 , function(i) {
    ee = translatedExps[[i]][ matchGff( allGff , translatedGffs[[i]] ) ]
    ee[is.na(ee)] = 0
    ee
  } )
  
  allCoverage = c(ct2607=14.7, ct553=7.5, sz581=13.2, sz641=15.6)
  expMatrix2 = expMatrix/allCoverage
  
  shared = apply( expMatrix2, 1, function(x) sum(x[3:4])>0 & sum(x[1:2])>0 )
  cased = apply( expMatrix2, 1, function(x) sum(x[3:4])>0 & sum(x[1:2])==0 )
  contd = apply( expMatrix2, 1, function(x) sum(x[1:2])>0 & sum(x[3:4])==0 )
  
  cols = rep(1,nrow(expMatrix))
  cols[shared] = 1
  cols[cased] = 2
  cols[contd] = 3
  
  de = t( apply( expMatrix2[shared,] , 1 , function(x) {
    cont = mean(x[1:2])
    case = mean(x[3:4])
    fc = case/cont
    p = t.test( x[1:2] , x[3:4] )$p.val
    c( case=case , cont=cont , fc=fc , p=p  )
  } ) )
  
  sum(de[,3]>2)
  sum(de[,3]<1/2)
  
  res = data.frame( total=nrow(expMatrix) , shared=sum(shared) , 
                    caseUnique=sum(cased) ,contUnique=sum(contd) , sharedCaseUpDE=sum(de[,3]>2) , 
                    sharedCaseDownDE=sum(de[,3]<1/2)  )
  
  outputDir = paste0(folder,"/comparisonResult","/",parai)
  system( paste( "mkdir -p" , outputDir ) )
  setwd(outputDir)
  
  write.table(res,"longReadDE.csv",col.names=T,row.names=F,quote=F,sep=",")

}




