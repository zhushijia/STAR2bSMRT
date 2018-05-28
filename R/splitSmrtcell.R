splitSmrtcell = function( alignments , outputDir , thres=0  )
{
  STAR2bSMRT.dir = system.file(package = "STAR2bSMRT")
  func = paste0("source " , STAR2bSMRT.dir , "/data/splitSmrtcell.sh")
  output = paste0(outputDir,"/smrtcells.uniq")
  sh1 = paste( func , "-uniq" , alignments , output )
  if( !file.exists(output) ) runSH(sh1)
  
  smrtcell = read.table(output,header=T)
  smrtcell = subset(smrtcell,count>=thres)
  for( i in seq_len(nrow(smrtcell)) )
  {
    system( paste0( "mkdir -p " , outputDir , "/separateSam/" , smrtcell[i,2] ) )
    outputI = paste0( outputDir, '/separateSam/' ,smrtcell[i,2] , "/Aligned.out.sam" )
    sh2 = paste0('grep "' , smrtcell[i,2] , '" ' , alignments ,' > ' , outputI )
    if( !file.exists(outputI) ) runSH(sh2)
  }
}

