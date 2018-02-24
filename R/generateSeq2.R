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
	print(genome)
	print(class(genome))
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
			sequence = list()
			for ( i in 1:nrow(exon) ) 
			  sequence[[i]] = substr( "adfasdfasdfadsfads" , 1 , 5 )
			
			cat("c\n")
			for ( i in 1:nrow(exon) ) 
			  sequence[[i]] = substr( as.character(genome[[ "chr2" ]]) , 1000000 , 1000010 )
			
			cat("d\n")
			for ( i in 1:nrow(exon) ) 
			  sequence[[i]] = subseq( genome[[ "chr2" ]] , exon$start[i] , exon$end[i] )
			
			cat("e\n")
			for ( i in 1:nrow(exon) ) 
			  sequence[[i]] = substr( genome[[ "chr2" ]] , exon$start[i] , exon$end[i] )
			
			cat("f\n")
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

