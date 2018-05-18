library(Biostrings)
ref = "/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa"
genome = readDNAStringSet(ref)
allSeq = lapply( gffs[1:4] , function(x) generateSeq( genome , isoform=x ) )

gffs2 = lapply( 1:4, function(i) {
  ind = which( allSeq[[i]]$translated )
  gffs[[i]][ind]
} )

exp2 = lapply( gffs2 , function(x) {
  sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
} )

names(gffs2) = names(gffs)[1:4]
names(exp2) = names(gffs)[1:4]


contGff2 =  unionGff(gffs2[[1]],gffs2[[2]])
caseGff2 =  unionGff(gffs2[[3]],gffs2[[4]])
allGff2 =  unionGff(caseGff2,contGff2)


allCoverage = c()
expMatrix = sapply( 1:4 , function(i) {
  ee = exp2[[i]][ matchGff( allGff2 , gffs2[[i]] ) ]
  ee[is.na(ee)] = 0
  ee
} )
expMatrix = expMatrix/allCoverage
colnames(expMatrix) = names(gffs)[1:4]



shared = apply( expMatrix, 1, function(x) all(x>0) )

de = t( apply( expMatrix , 1 , function(x) {
  cont = mean(x[1:2])
  case = mean(x[3:4])
  fc = case/cont
  p = t.test( x[1:2] , x[3:4] )$p.val
  c( case=case , cont=cont , fc=fc , p=p  )
} ) )


