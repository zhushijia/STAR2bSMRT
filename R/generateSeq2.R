#' generateSeq2
#'
#' generateSeq2 generates the transcript sequence
#' 
#' @param ref 
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
generateSeq2 = function( ref , isoform , exp , chrom , s , e )
{
	genome = readDNAStringSet(ref)
	cat("finish reading reference \n")
	dna = lapply( isoform , function(junc) {
	  cat("a\n")
	  chr = as.character(junc$chr[1])
	  junc$start = junc$start-1
	  junc$end = junc$end+1
	  junc = rbind( data.frame(chr=chr,start=0,end=s) , junc, data.frame(chr=chr,start=e,end=Inf) )
		exon = data.frame( chr=junc$chr[-1] , start=junc$end[-nrow(junc)] , end=junc$start[-1] )
		if( all((exon$end-exon$start)>0) )
		{
		  cat("b\n")
			sequence = sapply(1:nrow(exon), function(i) substr( genome[[ "chr2" ]] , exon$start[i] , exon$end[i] ) )
			cat("c\n")
			sequence = paste( sapply( sequence , function(x) as.character(x)) , collapse="")
			reverseComplement(DNAString(sequence))
		}
		
		} ) 
	
	dna = DNAStringSet( dna )
	names(dna) = paste0("SS",1:length(dna),"_exp",exp)
	# sum(as.numeric(sapply(fasta2,nchar)%%3)>0)

	suppressWarnings( protein <-translate(dna) )
	
	translated = sapply( protein , function(x) {
		y = as.character(x)
		grepl("^M",y) & grepl("[*]$",y) 
	} )
	
	list(dna=dna,protein=protein,translated=translated)
	
}

