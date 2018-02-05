phqvExp = function( phqv , outputDir )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/sh/phqvExp.sh")
	output = paste0(outputDir,"/phqv.exp")
	sh = paste( func , phqv , output  )
	if( !file.exists(output) ) runSH(sh)

	read.table( output , sep="\t" , header=T )
}
