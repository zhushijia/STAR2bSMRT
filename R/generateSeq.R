#' generateSeq
#'
#' @param genome DNAStringSet indicating the genome reference
#' @param isoform a list of data frame indicating the exon locations in isoforms
#'
#' @return a list comprsing the following items 
#' \itemize{
#'  \item {rna} character value indicating the RNA sequence of isoforms
#'  \item {protein} character value indicating the amino acid sequence of 
#'  isoforms
#'  \item {translated}  boolean value indicating whether the isoform can be 
#'  translated 
#' }
#' @export
#'
#' @examples
generateSeq = function( genome , isoform )
{
	
	rna = lapply( isoform , function(exon) {
	  
	  chr = as.character(exon$chr[1])
	  
	  if( all((exon$end-exon$start)>0) )
		{
			sequence = sapply(1:nrow(exon), function(i) subseq( genome[[ chr ]] , exon$start[i] , exon$end[i] ) )
			sequence = paste( sapply( sequence , function(x) as.character(x)) , collapse="")
			reverseComplement(DNAString(sequence))
		}
		
		} ) 
	
	rna = DNAStringSet( do.call(c,rna) )
	#names(rna) = names(isoform)
	# sum(as.numeric(sapply(fasta2,nchar)%%3)>0)

	suppressWarnings( protein <-translate(rna) )
	
	translated = sapply( protein , function(x) {
		y = as.character(x)
		grepl("^M",y) & grepl("[*]$",y) 
	} )
	
	list(rna=rna,protein=protein,translated=translated)
	
}

