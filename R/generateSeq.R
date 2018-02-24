#' generateSeq
#'
#' generateSeq generates the transcript sequence
#' 
#' @param genome 
#' @param isoform
#' @param exp
#' @param chrom 
#' @param s 
#' @param e 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
generateSeq = function( genome , isoform , exp , chrom , s , e , cores=1 )
{
	registerDoMC(cores)
  
	dna = list()
	for( i in 1:length(isoform) )
	{
	  cat(i,"\n")
	  junc = isoform[[i]]
	  cat(1,"\n")
	  chr = as.character(junc$chr[1])
	  cat(2,"\n")
	  junc$start = junc$start-1
	  junc$end = junc$end+1
	  cat(3,"\n")
	  junc = rbind( data.frame(chr=chr,start=0,end=s) , junc, data.frame(chr=chr,start=e,end=Inf) )
		
	  cat(4,"\n")
	  exon = data.frame( chr=junc$chr[-1] , start=junc$end[-nrow(junc)] , end=junc$start[-1] )
		if( all((exon$end-exon$start)>0) )
		{
		  cat(5,"\n")
		  sequence = sapply(1:nrow(exon), function(i) substr( genome[[ "chr2" ]] , exon$start[i] , exon$end[i] ) )
		  cat(6,"\n")
		  sequence = paste( sapply( sequence , function(x) as.character(x)) , collapse="")
		  cat(7,"\n")
		  dna[[i]] = reverseComplement(DNAString(sequence))
		}
	}
	
	cat(8,"\n")
  dna = DNAStringSet( dna )
  
	names(dna) = paste0("SS",1:length(dna),"_exp",exp)
	# sum(as.numeric(sapply(fasta2,nchar)%%3)>0)

	cat(9,"\n")
	suppressWarnings( protein <-translate(dna) )
	
	cat(10,"\n")
	translated = sapply( protein , function(x) {
		y = as.character(x)
		grepl("^M",y) & grepl("[*]$",y) 
	} )
	
	list(dna=dna,protein=protein,translated=translated)

	
}

