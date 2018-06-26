junc2tag = function(gff) 
{
	lapply( gff , function(read) 
	{
		tmp = paste( read$end[-nrow(read)]+1 , read$start[-1]-1 , sep="," )
		paste(tmp,collapse="; ")
	})
}

library(STAR2bSMRT,lib.loc="/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone4/setup")
library(foreach)
library(doMC)
registerDoMC(3)
LoutputDir = "/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/clone_2607_20171215/starLongNew"
LRalignment = paste0(LoutputDir,"/Aligned.out.sam")  
LRread = getReadByJI( LRalignment , LoutputDir )
LRinfo = getLRinfo( LRread )
cloneTag = lapply( LRinfo$LRtag[[1]] , function(x) paste( x[-length(x)], collapse="; " ) )

folder_2607="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/2607"
folder_result="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/result/STAR2bSMRT/pipelineNew_testParameters/comparisonResult"
parameters = dir(folder_2607)

for(i in 1:length(parameters))
{
	parai = paste0(folder_2607,"/",parameters[i])
	setwd(parai)
	
	gff = readGff(dir()[grepl(".gff$",dir())])
	tag = junc2tag(gff)
	matchedTopo = sapply( strsplit( names(gff)[ match(cloneTag,tag) ] ," "), function(x) gsub(";","",x[4]) )

	abundance = read.table( "output/abundance.tsv" , header=T )
	transcriptId = as.character(abundance$target_id)
	Lexp = sapply( strsplit(transcriptId,"_") , function(x) log(as.integer(gsub("exp","",x[2]))) )
	Sexp= log(abundance$tpm)
	translated = sapply( strsplit(transcriptId,"_") , function(x) x[3]=="translated" )
	info = data.frame( transcriptId , Lexp , Sexp=abundance$tpm , translated )

	col = rep(2,nrow(abundance))
	col[ which(transcriptId %in% matchedTopo) ] = 3
	col[ !translated ] = 1
	pch=rep(16,nrow(abundance))
	pch[!translated] = 17
	
	
	outputDir = paste0(folder_result,"/",parameters[i])
	system(paste0("mkdir -p ",outputDir))
	setwd(outputDir)
	pdf('topclone_2607.pdf')
	plot( Lexp , Sexp , pch=pch , col=col , main="Quantification by Long and Short Reads" ,  xlab="Log10 Long Read" , ylab="Log10 Short Read"  )
	abline(lm( Sexp~Lexp ))
	
	dev.off()

}

