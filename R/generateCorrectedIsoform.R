#' generateCorrectedIsoform
#' generates the corrected isoforms given the defined thresholds of ts and td
#' 
#' @param LRjunc a list of data frame, indicating the splicing junction sites 
#' obtained from long reads
#' @param SRjunc a list of data frame, indicating the splicing junction sites 
#' obtained from short reads
#' @param LRtag a list of character values, indicating the splicing junctions 
#' for each long read in the format of "junction_start,junction_end"
#' @param LRread a list of data frame, indicating the information of long read,
#' comprising id, chr, strand, start, end, junc
#' @param ts integer value indicating the number of short read support
#' @param td integer value indicating the distance between long and short read
#' junction sites
#' @param matchedLS the correspondence bewteen long read and short read-derived
#' junction sites
#' @param fuzzyMatch integer value indicating the tolerance distance  
#'
#' @return a list of corrected isoforms, comprising the following items:
#' \itemize{
#'  \item {isoform} a list of isoforms: data frame of junction sites
#'  \item {uniNum} integer value indicating the number of unique isoforms
#'  \item {frac} numeric value indicating the fraction of corrected reads 
#'  relative to total read counts
#'  \item {isoformCount} integer value indicating the read count for each isoform
#'  \item {LSjuncCount} data frame representing the information of each isoform,
#'  comprising columns: srCount, chr, start, end, motif, and lrCount
#'  \item {normalizedIsoformCount} numeric values indicating the normalized 
#'  read count for each isoform when integrating multiple smrtcells
#' }
#' 
#' 
#' @export
#'
#' @examples
#' 
#' 
generateCorrectedIsoform = function( LRjunc , SRjunc, LRtag , LRread  , ts , td , matchedLS=NULL , fuzzyMatch=100 )
{
  
  CHR = intersect( names(SRjunc) , names(LRjunc) )
  
  allCorrectIsoform = foreach( k=1:length(CHR) ) %dopar% {
    
    chr = CHR[k]
    lrc = LRjunc[[chr]]
    src = SRjunc[[chr]]
    tag = LRtag[[chr]]
    read = LRread[[chr]]
    
    if( !is.null(matchedLS) )
    {
      SRmatch = matchedLS[[chr]]
      index = which( SRmatch[,'LSdis']<=td & src$count[SRmatch[,'SRindex']]>=ts ) 
    } else {
      src = subset( src , count>=ts )
      SRmatch = matchLSjuncOneChr( lrc , src ,fuzzyMatch )
      index = which( SRmatch[,'LSdis']<=td ) 
    }
    
    #LRcorres = data.frame(lrc[index,],SRmatch[index,])
    LRcorres = data.frame(SRmatch,lrc)[index, ]
    correctTag = paste(LRcorres$start , LRcorres$end , sep="," )
    range = sapply( tag , function(x) all(x%in%correctTag) )
    isoform = sapply( tag[range] , function(x) paste( sort(subset(LRcorres,correctTag%in%x)$SRindex) , collapse="_" ) )
    
    #correctIsoformExp = tapply( read$coverage[range] , isoform , sum )
    #srindex = lapply( strsplit( names(correctIsoformExp) , "_" ) , function(x) as.integer(x) )
    
    info = data.frame( read[range,] , isoform  )
    correctIsoformCount = cast( info , isoform ~ group, sum , value="coverage")
    srindex = lapply( strsplit( as.character(correctIsoformCount$isoform) , "_" ) , 
                      function(x) as.integer(x) )
    
    correctIsoform = lapply( srindex , function(ind) { 
      x=src[ ind , c('chr','start','end') ]
      x[ order(x$start), ]
    } )
    
    
    #lrCount = tapply( LRcorres$count , LRcorres$SRindex , sum )
    #LSjuncCount = data.frame( lrCount , src[as.integer(names(lrCount)),] )
    #colnames(LSjuncCount) = c("lrCount","srCount","chr","start","end","motif")
    
    lrCount = apply( as.matrix(LRcorres[,-c(1:5) ]) , 2 , 
                     function(x) tapply( x , LRcorres$SRindex , sum ) )
    LSjuncCount = data.frame( src[as.integer(rownames(lrCount)),], lrCount )
    colnames(LSjuncCount)[1] = c("srCount")
    
    correctIsoformNum = length(unique(isoform))
    
    totalReadCount = tapply(read$coverage, read$group, sum)
    groupReadCount = apply(as.matrix(correctIsoformCount[,-1]),2,sum)
    correctIsoformFrac = groupReadCount/totalReadCount
     
    list( isoform=correctIsoform, uniNum=correctIsoformNum, 
          frac=correctIsoformFrac, isoformCount=correctIsoformCount, 
          LSjuncCount=LSjuncCount )
    
  }
  
  names(allCorrectIsoform) = CHR
  
  LSjuncCount = do.call( rbind , 
                            lapply( allCorrectIsoform, function(x) x$LSjuncCount ) )
  
  if(ncol(LSjuncCount)>6)
  {
    y = log(LSjuncCount[,1]+1)
    x = as.matrix( log(LSjuncCount[,-c(1:5)]+1) )
    posCoef = as.matrix( coef(nnnpls(x,y,con=rep(1,ncol(x)))) )
  } else {
    posCoef = as.matrix(1)
  }
  
  for( chr in CHR )
  {
    x = as.matrix(allCorrectIsoform[[chr]][['isoformCount']][,-1])
    allCorrectIsoform[[chr]][['normalizedIsoformCount']] = as.numeric( x %*% posCoef)
  }
  
  allCorrectIsoform
  
}
