#' generateSeq
#'
#' @param genome 
#' @param isoform 
#'
#' @return
#' @export
#'
#' @examples
generateSeq = function( genome , isoform )
{
	
	dna = lapply( isoform , function(exon) {
	  
	  chr = as.character(exon$chr[1])
	  
	  if( all((exon$end-exon$start)>0) )
		{
			sequence = sapply(1:nrow(exon), function(i) subseq( genome[[ chr ]] , exon$start[i] , exon$end[i] ) )
			sequence = paste( sapply( sequence , function(x) as.character(x)) , collapse="")
			reverseComplement(DNAString(sequence))
		}
		
		} ) 
	
	dna = DNAStringSet( do.call(c,dna) )
	#names(dna) = names(isoform)
	# sum(as.numeric(sapply(fasta2,nchar)%%3)>0)

	suppressWarnings( protein <-translate(dna) )
	
	translated = sapply( protein , function(x) {
		y = as.character(x)
		grepl("^M",y) & grepl("[*]$",y) 
	} )
	
	list(dna=dna,protein=protein,translated=translated)
	
}

