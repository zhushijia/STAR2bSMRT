starShort = function( genomeDir , SR1 , SR2 , outputDir )
{
	STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
	func = paste0("source " , STAR2bSMRT.dir , "/sh/starShort.sh")
	sh = paste( func , genomeDir , SR1 , SR2 , outputDir  )
	runSH(sh)
}
