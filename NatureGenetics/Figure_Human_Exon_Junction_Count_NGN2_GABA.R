
calculateExonExp = function(cc_fracs,cc_exps)
{
  cc_exps = as.matrix(cc_exps)
  exonExpExpEachTranscript = apply( cc_fracs , 2 , function(x) x*rowSums(cc_exps) )
  colSums(exonExpExpEachTranscript)
}

calculateJuncFrac = function(cc_fracs)
{
  isoform = apply(cc_fracs , 1 , function(x) {
    iso1 = paste( which(x>0), collapse = "_" )
    paste0( "_" , iso1 , "_")
  })
  
  junc = matrix(ncol=ncol(cc_fracs),nrow=ncol(cc_fracs),data=0)
  colnames(junc) = colnames(cc_fracs)
  rownames(junc) = colnames(cc_fracs)
  
  for(i in 1:ncol(cc_fracs))
  {
    for(j in i:ncol(cc_fracs))
    {
      if(i==j)
      {
        junc[i,i] = sum( grepl( paste("",i,"",sep="_") , isoform ) )
      } else {
        junc[i,j] = sum( grepl( paste("",i,j,sep="_") , isoform ) )
      }
    }
    
  }
  
  frac = junc/diag(junc)
  frac[lower.tri(frac)] = -1
  diag(frac) = 0
  frac[is.na(frac)] = 0
  frac
  
}

##########################################################################################
#################   load codes for heatmap and gff
##########################################################################################

source("SourceCode_Gff_Op.R")
source("SourceCode_Heatmap.R")

##########################################################################################
#################   human samples and basic information
##########################################################################################
samples= c("553_NGN2", "553_GABA") 
human_annotation = read.table('Data/human_NRXN1_Exon_Annotation_hg19.txt',sep='\t',header=T)
human_nrxn1_chr = "chr2"
human_nrxn1_start = 50149082 
human_nrxn1_end = 51255411
groups = as.list(1:length(samples))
names(groups) = samples
threshold=7

# filter translated isoforms
gffs = lapply( samples , function(sample)
{
    gff = paste0("Data/Human_NRXN1alpha_Translated_Thres7_",sample,".gff")
    readGff( gff , chrom=human_nrxn1_chr , s=human_nrxn1_start , e=human_nrxn1_end )
} )
names(gffs) = samples


translated = lapply( gffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) gsub("exp|;","",p[5])=="translated" )
}   )

translatedGffs = lapply( 1:length(gffs) , function(i) gffs[[i]][ translated[[i]] ]  )
names(translatedGffs) = names(gffs)

translatedExps = lapply( translatedGffs , function(x) {
    sapply( strsplit(names(x),'_'), function(p) as.integer(gsub("exp|;","",p[4])) )
}   )
names(translatedExps) = names(gffs)

translatedGffs = lapply( 1:length(translatedGffs) , function(i) translatedGffs[[i]][ translatedExps[[i]]>=threshold ] )
translatedExps = lapply( 1:length(translatedGffs) , function(i) translatedExps[[i]][ translatedExps[[i]]>=threshold ] )
names(translatedGffs) = names(gffs)
names(translatedExps) = names(gffs)

translatedFracs = lapply( translatedGffs, function(gff) {
    frac = getFrac(gff, human_annotation)
    as.data.frame(round(frac))
} )

translatedExonExp = lapply( 1:length(gffs) , function(i) 
    calculateExonExp(translatedFracs[[i]],translatedExps[[i]]))
names(translatedExonExp) = names(gffs)

translatedJuncFrac = lapply( 1:length(gffs) , function(i) 
    calculateJuncFrac(translatedFracs[[i]]) )
names(translatedJuncFrac) = names(gffs)

##########################################################################################
############# draw figure for 553_NGN2
##########################################################################################
dir.create("Figs")

pdf("Figs/Fig3B_Human_Exon_Junction_Count_NGN2_GABA.pdf")

i = 1
juncFrac = translatedJuncFrac[[i]]
juncExp = translatedExonExp[[i]]
Name = reorder( rownames(juncFrac) , 1:nrow(juncFrac) )
X = data.frame(Name,juncFrac)
Y = data.frame(Name,y=juncExp)
heatmapGradientColor2(X,Y,names(translatedJuncFrac)[i])


##########################################################################################
############# draw figure for 553_GABA
##########################################################################################
i = 2
juncFrac = translatedJuncFrac[[i]]
juncExp = translatedExonExp[[i]]
Name = reorder( rownames(juncFrac) , 1:nrow(juncFrac) )
X = data.frame(Name,juncFrac)
Y = data.frame(Name,y=juncExp)
heatmapGradientColor2(X,Y,names(translatedJuncFrac)[i])

dev.off()


