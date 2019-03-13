getPSI = function(ef)
{
  library(STAR2bSMRT,lib.loc="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/code/STAR2bSMRT/githubClone5/setup")
  annotation = read.table('/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/ToolCompare/NRXN1_hg19_ExonAnnotations_shijia.txt',sep='\t',header=T)
  parent="/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping"
  chrom = "chr2"
  s = 50147488
  e = 51259537
  SoutputDir = paste0(parent,'/',ef)
  SRalignment = paste0(SoutputDir,"/alignments.bam")
  #SRjunc = getJuncBySam( SRalignment, SoutputDir, chrom=chrom, s=s, e=e )
  SRjunc = getJuncBySJout( SJout="SJ.out.tab",outputDir=SoutputDir,chrom=chrom,s=s,e=e )
  junc = SRjunc[['chr2']]
  
  jstop = function(exon) { subset(annotation , Exon==exon)$fullStop-1 }
  jstart = function(exon) { subset(annotation , Exon==exon)$fullStart+1 }
  PSI = function( exonLeft , exonMid , exonRight )
  {
    LeftIn = with( junc, count[ end == jstop(exonLeft) & start == jstart(exonMid) ] )
    RightIn = with( junc, count[ end == jstop(exonMid) & start == jstart(exonRight) ] )
    Ex = with( junc, count[ end == jstop(exonLeft) & start == jstart(exonRight) ] )
    
    if(length(LeftIn)==0) LeftIn=0
    if(length(RightIn)==0) RightIn=0
    if(length(Ex)==0) Ex=0
    if (exonLeft=='exon22')
      Ex = Ex + with( junc, count[ end == jstop('exon20') & start == jstart(exonRight) ] )
    
    psi = (LeftIn+RightIn)/(LeftIn+RightIn+Ex)
    psi
  }
  
  psi = c()
  psi['3a'] = PSI('exon2','exon3a','exon4')
  psi['3b'] = PSI('exon2','exon3b','exon4')
  psi['4'] = PSI('exon3a','exon4','exon5')
  psi['5'] = PSI('exon4','exon5','exon6')
  
  psi['7a'] = PSI('exon6','exon7a','exon8')
  psi['7b'] = PSI('exon6','exon7b','exon8')
  
  psi['12'] = PSI('exon11','exon12','exon13')
  psi['17'] = PSI('exon16','exon17','exon18')
  psi['21'] = PSI('exon20','exon21','exon22')
  psi['23a'] = PSI('exon22','exon23a','exon24')
  psi['23b'] = PSI('exon22','exon23b','exon24')
  
  psi
  
}


setwd("/sc/orga/projects/schzrnas/sjzhu/Project/NRXN/data/TargetShortReadRerun_Genewiz/MH1804035/mapping")

PBS_2607 = getPSI('EF15')
KCL_2607 = getPSI('EF17')
PBS_553 = getPSI('EF21')
KCL_553 = getPSI('EF23')
PBS_641 = getPSI('EF18')
KCL_641 = getPSI('EF20')

KCL_case = KCL_641
PBS_case = PBS_641
ratio_cont_2607 = KCL_2607/PBS_2607
ratio_cont_553 = KCL_553/PBS_553
ratio_cont = cbind( ratio_cont_2607 , ratio_cont_553 )

#KCL_cont = ( KCL_2607 + KCL_553 ) / 2
#PBS_cont = ( PBS_2607 + PBS_553 ) / 2

ratio_case_mean = KCL_case/PBS_case
ratio_case_se = rep(0,length(ratio_case_mean))
ratio_cont_mean = apply(ratio_cont,1,mean)
ratio_cont_se = apply(ratio_cont,1,function(x) sd(x)/sqrt(2)  )

PSI_ratio = c( matrix( rbind(ratio_cont_mean,ratio_case_mean,0) , nrow=1 ))
PSI_error = c( matrix( rbind(ratio_cont_se,ratio_case_se,0) , nrow=1 ))
names(PSI_ratio) = c( matrix(rbind(names(ratio_case_mean),"",""), nrow=1 ))

pdf("PSI_FC_KCLvsPSB_byJunc_withError.pdf")
barx <- barplot(PSI_ratio,col=rep(c('grey','darkred','black'),11),las=3,ylab='PSI fold change',ylim=c(0,1.4))
arrows(barx , PSI_ratio+PSI_error, barx , PSI_ratio, angle=90, code=3, length=0.05)
dev.off()
