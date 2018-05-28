getReadByCigar = function( alignments , outputDir )
{
  
  thres = 20
  sam <- file( alignments, "r", blocking = FALSE)
  
  read = readLines(sam,n=1)
  while( length(read)>0 )
  {
    cat(read,"\n")
    sep = strsplit(read,"\t")[[1]]
    if(length(sep)<10)
    {
      read = readLines(sam,n=1)
    } else {
      break
    }
  }
  
  juncs = list()
  
  id = c()
  chr = c()
  strand = c()
  start = c()
  end = c()
  junc = c()
  
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
      index = which( cigar[,1]>thres & as.character(cigar[,2]) %in% c("N","D") )
      
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
        
        id[num] = sep[1]
        chr[num] = sep[3]
        strand[num] = NA
        start[num] = as.integer(sep[4]) - 1 # to match bamToBed
        end[num] = start[num] + regionLen 
        junc[num] = paste(apply(data.frame(juncL,juncR),1,function(x)paste(x,collapse=",")) , collapse="," )
  
      }
      
      
    }
    
    read = readLines(sam,n=1)
    
  }
  
  close(sam)
  
  alignments.read = data.frame( id , chr , strand , start , end , junc )
  
  write.table( alignments.read , paste0(outputDir,"/alignments.read") , col.names=TRUE , row.names=FALSE , sep="\t" , quote=FALSE )
  
  return(alignments.read)
  
}



