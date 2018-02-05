getJunc = function( alignments , outputDir , chrom=NULL ,  s , e )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/sh/getJunc.sh")
	output = paste0(outputDir,"/alignments.junc")
	sh = paste( func , alignments , output  )
	if( !file.exists(output) ) runSH(sh)

	junc = read.table( output , sep="\t" , header=T )

	if( !is.null(chrom) )
	{
		junc = subset( junc , chr %in% chrom & start>=s & end<=e )
	}
	
	split( junc , as.character(junc$chr) )

}
