##########################################################################################
#################   load codes for heatmap and gff
##########################################################################################

source("SourceCode_Gff_Op.R")
source("SourceCode_Heatmap.R")

##########################################################################################
#################   load codes for heatmap and gff
##########################################################################################
samples = c("p3Del1","p3Del2","Cont1","Cont2","Adult1","Adult2","Adult3")
case_ind = c(1,2)
cont_ind = c(3,4)
valid_ind = c(5:7)
human_annotation = read.table('Data/human_NRXN1_Exon_Annotation_hg19.txt',sep='\t',header=T)
case_annotation = human_annotation
cont_annotation = human_annotation
human_nrxn1_chr = "chr2"
human_nrxn1_start = 50149082 
human_nrxn1_end = 51255411
threshold=7

##########################################################################################
#################   load codes for heatmap and gff
##########################################################################################

# filter translated isoforms
gffs = lapply( samples , function(sample)
{
    gff = paste0("Data/Human_NRXN1alpha_Translated_Thres7_",sample,".gff")
    readGff( gff , chrom=human_nrxn1_chr , s=human_nrxn1_start , e=human_nrxn1_end )
} )
names(gffs) = samples

translated = lapply( gffs , function(x) {
sapply( strsplit(names(x),'_'), function(p) p[5]=="translated;" )
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

# union and intersect samples 
caseGff = NULL
for (i in case_ind) caseGff = unionGff(caseGff,translatedGffs[[i]])
contGff = NULL
for (i in cont_ind) contGff = unionGff(contGff,translatedGffs[[i]])

case_fracs = getFrac(caseGff,case_annotation)
cont_fracs = getFrac(contGff,cont_annotation)

case_fracs2 = as.data.frame(round(case_fracs))
cont_fracs2 = as.data.frame(round(cont_fracs))

range = matchGff(contGff,caseGff)
range = range[!is.na(range)]

caseUniFrac = -case_fracs2[-range, ] # -1 means case unique 
commonFrac = case_fracs2[range, ]

caseUniGff = caseGff[-range]
commonGff = caseGff[range]

range = matchGff(caseGff,contGff)
range = range[!is.na(range)]

contUniFrac = 2*cont_fracs2[-range, ]
contUniGff = contGff[-range]

unique(do.call(c,caseUniFrac))
unique(do.call(c,contUniFrac))
unique(do.call(c,commonFrac))

fracs = rbind(caseUniFrac,contUniFrac,commonFrac)

allGff = c( caseUniGff, contUniGff , commonGff )
exps = sapply( 1:length(translatedGffs) , function(i) {
  ee = translatedExps[[i]][ matchGff( allGff , translatedGffs[[i]] ) ]
  ee[is.na(ee)] = 0
  ee
} )

# normalized to sum
normalizedFactor = colSums(exps)/min(colSums(exps[,c(case_ind,cont_ind)]))
exps = sapply( 1:ncol(exps) , function(i) exps[,i]/normalizedFactor[i] )
colnames(exps) = names(translatedGffs)
exps = data.frame(exps)
de = as.numeric( apply( exps , 1 , function(x) {
    case = mean(x[case_ind])
    cont = mean(x[cont_ind])
    log2(case+1)-log2(cont+1)
} ) )


Name = reorder( paste0('isoform',1:nrow(exps)) , de )
logExps = data.frame( Name, log10(exps[,c(case_ind,cont_ind)]+1) )
fracs = data.frame( Name , fracs ) 

################################################################################
########## validate by adult sample
################################################################################

validGff = NULL
for (i in valid_ind) validGff = unionGff(validGff,translatedGffs[[i]])
validate = rep(0,length(allGff))
validate[ matchGff(validGff,allGff) ] = 1
validate2 = data.frame( Name , validate ) 
validateName = paste( c("validateBy", names(gffs)[valid_ind]), collapse="_")

################################################################################
#############  isoform schematics
################################################################################
pdf("Figs/Fig4A_Human_Isoform_Schematics_hiPSC_case_cont.pdf")
# draw isoform schematics
heatmapMannualColor3(fracs, case_col="darkred", cont_col="darkgrey", shared_col="darkorange")
# draw expression
heatmapGradientColor1(logExps)
# draw validated isoforms
heatmapMannualValidateColor(validate2)
dev.off()


