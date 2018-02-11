#' generateSeq
#'
#' @param fasta 
#' @param Junc 
#' @param chrom 
#' @param s 
#' @param e 
#'
#' @return
#' @export
#'
#' @examples
generateSeq = function( fasta , Junc , chrom , s , e )
{
	library(Biostrings)
	#genome = readDNAStringSet( fasta )
	genome = readDNAStringSet("/hpc/users/zhus02/schzrnas/sjzhu/RNAseq/Reference/hg19/reference/hg19.fa")


	fasta = lapply( correctJunc , function(s) {
		s$start = s$start-1
		s$end = s$end+1
		s = rbind( data.frame(chr="chr2",start=0,end=50149082) , s, data.frame(chr="chr2",start=51255411,end=Inf) )
		exon = data.frame( chr=s$chr[-1] , start=s$end[-nrow(s)] , end=s$start[-1] )
		if( all((exon$end-exon$start)>0) )
		{
			seq = sapply(1:nrow(exon), function(i) substr( genome[[ "chr2" ]] , exon$start[i] , exon$end[i] ) )
			paste( sapply( seq,function(x) as.character(x)) , collapse="")
		}
		
		} )
	names(fasta) = paste0("seq",1:length(fasta),"_exp",correctExp)
	sort(sapply(fasta,nchar))
	fasta = fasta[ !is.na(as.integer(sapply(fasta,nchar))) ]
	fasta2 = DNAStringSet( lapply(fasta,function(x) reverseComplement(DNAString(x)) ) )
	names(fasta2) = paste0("seq",1:length(fasta),"_exp",correctExp)

	as.numeric(sapply(fasta2,nchar)/3)
	fa=paste0("test_i",i,"_j",j,"_2")
	folder=paste0("/hpc/users/zhus02/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/LongReads/Smrtportal_24463_2607/starLongNewUncorrectedHighQVConsensus/STAR2bSMRT/",fa)
	dir.create(folder)
	setwd(folder)
	#writeXStringSet( fasta2,paste0(fa,".fa"))
	translate(fasta2[[1]])

	translation = translate(fasta2)
	se = sapply( translation , function(x) {
		y = as.character(x)
		grepl("^M",y) & grepl("[*]$",y) 
	} )
	
	num = sapply( translation , function(x) {
		y = as.character(x)
		length(gregexpr("[*]",y)[[1]])
	} )
	
	sum(num>=1)
	sum(se & num==1)
	length(translation)
	
	sum(correctExp[ se & num==1 ])
	sum(correctExp)
	
}

