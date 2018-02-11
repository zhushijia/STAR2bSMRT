getRead = function( alignments , outputDir )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/data/getRead.sh")
	output = paste0(outputDir,"/alignments.read")
	sh = paste( func , alignments , output  )
	if( !file.exists(output) ) runSH(sh)

	read.table( output , sep="\t" , header=T )
}

