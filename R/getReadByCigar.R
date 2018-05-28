#' getReadByCigar
#'
#' @param alignments 
#' @param outputDir 
#' @param intronMinSize 
#'
#' @return
#' @export
#'
#' @examples
#' 
getReadByCigar = function( alignments , outputDir , intronMinSize=20 )
{
  output <- paste0( outputDir , "/alignments.read" )
  
  if( !file.exists(output) )
  {
    samConn <- file( alignments, "r", blocking = FALSE)
    
    outputTmp <- paste0(outputDir,"/tmp.read")
    write( paste("id","chr","strand","start","end","junc",sep="\t") , outputTmp )
    
    read = readLines(samConn,n=1)
    while( length(read)>0 )
    {
      cat(read,"\n")
      sep = strsplit(read,"\t")[[1]]
      if(length(sep)<10)
      {
        read = readLines(samConn,n=1)
      } else {
        break
      }
    }
    
    num=0
    
    while( length(read)>0 )
    {
      sep = strsplit(read,"\t")[[1]]
      
      if( nchar(sep[6])>1 )
      {
        s = as.integer(sep[4])
        cigar = parseCigar(sep[6])
        #readLen = with( cigar , sum( num[ Op%in%c('M','I','S','=','X') ] ) )
        cigar = subset( cigar , !Op%in%c("I","S") )
        regionLen = sum( cigar$num )
        index = which( cigar[,1]>intronMinSize & as.character(cigar[,2]) %in% c("N","D") )
        
        if( length(index) > 0 )
        {
          num = num + 1
          cat(num,"\n")
          juncL = c()
          juncR = c()
          
          for(i in 1:length(index))
          {
            juncL[i] = s + sum(cigar$num[1:(index[i]-1)])
            juncR[i] = s + sum(cigar$num[1:index[i]]) - 1
          }
          
          id = sep[1]
          chr = sep[3]
          strand = NA
          start = as.integer(sep[4]) - 1 # to match bamToBed
          end = start + regionLen 
          junc = paste( apply( data.frame(juncL,juncR), 1,
                              function(x) paste( x, collapse=",")) , collapse=",")
          
          res = paste( id , chr , strand , start , end , junc , sep="\t" )
          write( res , outputTmp , append=TRUE )
        }
      }
      
      read = readLines(samConn,n=1)
      
    }
    
    close(samConn)
    file.rename(outputTmp , output )
  }
  
  alignments.read <- read.table( output , sep="\t" , header=TRUE )
  return(alignments.read)
  
}



