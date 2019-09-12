library(VennDiagram)


readGff <- function( file , chrom=NULL , s=0 , e=Inf )
{
    
    gff = read.table(file,sep='\t')
    colnames(gff)[c(1,3:5,7)] = c('chr','type','start','end','strand')
    exon = subset( gff , type=='exon' & start>=s & end<=e )
    
    if( !is.null(chrom) )
        exon = subset( exon , chr==chrom )
    
    split( exon[,c(1,4:5)] , exon[,9] )
    
}

writeGff <- function( isoform , file = "" )
{
    gff = list()
    for (i in 1:length(isoform))
    {
        exoni = isoform[[i]]
        len = nrow(exoni)
        chr = rep( as.character(exoni$chr[1]) , len+1 )
        generator = rep("SS",len+1)
        type = c( 'transcript' , rep( 'exon',len ) )
        start = c( exoni$start[1] , exoni$start )
        end = c( exoni$end[len] , exoni$end )
        strand = "-"
        id = paste0( 'gene_id "SS.1"; transcript_id "',names(isoform)[i],'";' )
        gff[[i]] = data.frame( chr,generator,type,start,end,".",strand,'.',id )
    }
    
    gff = do.call(rbind,gff)
    
    write.table( gff , file , col.names=F , row.names=F , sep="\t" , quote=F )
    
}

matchGff <- function( gff1 , gff2 )
{
    g1 = lapply( gff1,function(x) paste( x$chr[1] , paste( paste( x$start , x$end  ),collapse="; "),sep=": " ) )
    g2 = lapply( gff2,function(x) paste( x$chr[1] , paste( paste( x$start , x$end  ),collapse="; "),sep=": " ) )
    match(g1,g2)
}

unionGff <- function( gff1 , gff2 )
{
    index = matchGff(gff1 , gff2)
    index = index[!is.na(index)]
    index2 = setdiff( seq_along(gff2) , index  )
    c( gff1 , gff2[index2] ) 
}

vennDiagramGff <- function(gffList,filename)
{
    toTag = function(x) paste( x$chr[1] , paste( paste( x$start , x$end  ),collapse="; "),sep=": " )
    tagList = lapply( gffList , function(y) sapply( y , toTag ) )
    tagAll = unique( do.call(c, tagList ))
    indexList <- lapply( tagList , function(x) which(tagAll%in%x) )
    names(indexList) = names(gffList)
    venn.diagram(indexList, filename=filename , height = 2000, width = 2000, fill=1:length(gffList) , alpha=0.4 ,imagetype = "png")
}

gffInfo = function(gff, annotation, threshold)
{
    translated = sapply( strsplit(names(gff),'_'), function(p) grepl("^translated",p[5]) )
    translatedGff = gff[ translated ] 
    translatedExp = sapply( strsplit(names(translatedGff),'_'), function(p) as.integer(gsub("exp","",p[4])) )
    
    translatedGff = translatedGff[ translatedExp>=threshold ]
    translatedExp = translatedExp[ translatedExp>=threshold ]
    
    translatedFrac = getFrac(translatedGff, annotation)
    translatedFrac2 = as.data.frame(round(translatedFrac))
    translatedTag = apply(translatedFrac2,1,function(x) paste(colnames(translatedFrac2)[x>0],collapse="_") )
    list(gff=translatedGff, exp=translatedExp, tag=translatedTag,frac=translatedFrac2)
}


getFrac = function(gff,annotation)
{
    
    frac = function(gs,annot_sites)
    {
        sapply(annot_sites,function(x) {
            max( sapply(gs,function(y) mean(x%in%y)) )
        } )
    }
    
    gff_sites = lapply(gff,function(y) apply(y,1,function(z) z[2]:z[3] ) )
    annot_sites = apply(annotation,1,function(z) z[2]:z[3] )
    fracs = do.call(rbind,lapply( gff_sites , function(gs) frac(gs,annot_sites) ) )
    colnames(fracs) = as.character(annotation$Exon)
    rownames(fracs) = NULL
    fracs = data.frame(fracs)
    
    fracs$exon3b[ fracs$exon3a==1 & fracs$exon3b==1 ] = 0
    fracs$exon7a[ fracs$exon7a==1 & fracs$exon7b==1 ] = 0
    fracs$exon23a[ fracs$exon23a==1 & fracs$exon23b==1 ] = 0
    
    fracs
    
}
