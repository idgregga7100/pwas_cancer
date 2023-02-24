#false discovery rate instead of bonferroni (less stringent), test within each phenotype instead of all 7 combined
#p.adjust(p, method=p.adjust.methods, n=length(p))
#methods options are 'holm' 'hochberg' 'hommel' 'bonferroni' 'BH' 'BY' 'fdr' 'none' we want 'fdr' this time
#/home/isabelle/topmed/proteome/Schlidkraut/SPred_out/*_PCAIR_baseline_models_rho0.1_zpval0.05_*.csv ALL/CAU, 7 phenos
#pvalue is pvalue column
#/home/isabelle/topmed/proteome/Schlidkraut/INTERVAL_predixcan/INTERVAL_PrediXcan_*.csv
#pvalue is pvalue column
#/home/isabelle/FUSION/fusion_twas-master/cancer_pwas/FUSION_PWAS_*.csv
#TWAS.P is pvalue column

library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

phenos<-c('breast','breast_erpos','breast_erneg','endo','ovary','ovary_hgsoc','prostate')
pops<-c('ALL','CAU')

for (pop in pops){
  for (pheno in phenos){
    pwas<-fread('/home/isabelle/topmed/proteome/Schlidkraut/SPred_out/'%&%pop%&%'_PCAIR_baseline_models_rho0.1_zpval0.05_'%&%pheno%&%'.csv',stringsAsFactors = F)
    fdr_pvalue<-p.adjust(pwas$pvalue, method = c('fdr'))
    pwasfdr<-cbind(pwas,fdr_pvalue)
    fwrite(pwasfdr, '/home/isabelle/topmed/proteome/Schlidkraut/'%&%pop%&%'_baseline_rhozpval_'%&%pheno%&%'_FDR.txt',sep = '\t',quote = F,row.names = F)
  }
}

for (pheno in phenos){
  pwas<-fread('/home/isabelle/topmed/proteome/Schlidkraut/INTERVAL_predixcan/INTERVAL_PrediXcan_'%&%pheno%&%'.csv',stringsAsFactors = F)
  fdr_pvalue<-p.adjust(pwas$pvalue, method = c('fdr'))
  pwasfdr<-cbind(pwas,fdr_pvalue)
  fwrite(pwasfdr, '/home/isabelle/topmed/proteome/Schlidkraut/INTERVAL_PrediXcan_'%&%pheno%&%'_FDR.txt', sep = '\t', quote = F,row.names = F)
}

for (pheno in phenos){
  pwas<-fread('/home/isabelle/FUSION/fusion_twas-master/cancer_pwas/FUSION_PWAS_'%&%pheno%&%'.csv',stringsAsFactors = F)
  fdr_pvalue<-p.adjust(pwas$TWAS.P, method = c('fdr'))
  pwasfdr<-cbind(pwas,fdr_pvalue)
  fwrite(pwasfdr, '/home/isabelle/FUSION/FUSION_PWAS_'%&%pheno%&%'_FDR.txt', sep = '\t', quote = F,row.names = F)
}

#do a big combined pheno file for ease of joining later

combopwas<-fread('/home/isabelle/topmed/proteome/Schlidkraut/CAU_baseline_rhozpval_breast_FDR.txt',stringsAsFactors=F)
Phenotype<-rep_len(c('breast'),length.out=(nrow(combopwas)))%>%as.data.frame()
combopwas$Phenotype<-Phenotype

phenos<-c('breast_erpos','breast_erneg','endo','ovary','ovary_hgsoc','prostate')
for (pheno in phenos){
  pwas<-fread('/home/isabelle/topmed/proteome/Schlidkraut/CAU_baseline_rhozpval_'%&%pheno%&%'_FDR.txt',stringsAsFactors=F)
  Phenotype<-rep_len(pheno,length.out=(nrow(pwas)))%>%as.data.frame()
  pwas$Phenotype<-Phenotype
  combopwas<-rbind(combopwas,pwas)
}
fwrite(combopwas,'/home/isabelle/topmed/proteome/Schlidkraut/CAU_baseline_rhozpval_combined_FDR.txt',sep='\t',quote=F,row.names=F)

combopwas<-fread('/home/isabelle/topmed/proteome/Schlidkraut/ALL_baseline_rhozpval_breast_FDR.txt',stringsAsFactors=F)
Phenotype<-rep_len(c('breast'),length.out=(nrow(combopwas)))%>%as.data.frame()
combopwas$Phenotype<-Phenotype

phenos<-c('breast_erpos','breast_erneg','endo','ovary','ovary_hgsoc','prostate')
for (pheno in phenos){
  pwas<-fread('/home/isabelle/topmed/proteome/Schlidkraut/ALL_baseline_rhozpval_'%&%pheno%&%'_FDR.txt',stringsAsFactors=F)
  Phenotype<-rep_len(pheno,length.out=(nrow(pwas)))%>%as.data.frame()
  pwas$Phenotype<-Phenotype
  combopwas<-rbind(combopwas,pwas)
}
fwrite(combopwas,'/home/isabelle/topmed/proteome/Schlidkraut/ALL_baseline_rhozpval_combined_FDR.txt',sep='\t',quote=F,row.names=F)

combopwas<-fread('/home/isabelle/topmed/proteome/Schlidkraut/INTERVAL_PrediXcan_breast_FDR.txt',stringsAsFactors=F)
Phenotype<-rep_len(c('breast'),length.out=(nrow(combopwas)))%>%as.data.frame()
combopwas$Phenotype<-Phenotype

phenos<-c('breast_erpos','breast_erneg','endo','ovary','ovary_hgsoc','prostate')
for (pheno in phenos){
  pwas<-fread('/home/isabelle/topmed/proteome/Schlidkraut/INTERVAL_PrediXcan_'%&%pheno%&%'_FDR.txt',stringsAsFactors=F)
  Phenotype<-rep_len(pheno,length.out=(nrow(pwas)))%>%as.data.frame()
  pwas$Phenotype<-Phenotype
  combopwas<-rbind(combopwas,pwas)
}
fwrite(combopwas,'/home/isabelle/topmed/proteome/Schlidkraut/INTERVAL_PrediXcan_combined_FDR.txt',sep='\t',quote=F,row.names=F)

combopwas<-fread('/home/isabelle/FUSION/FUSION_PWAS_breast_FDR.txt',stringsAsFactors=F)
Phenotype<-rep_len(c('breast'),length.out=(nrow(combopwas)))%>%as.data.frame()
combopwas$Phenotype<-Phenotype

phenos<-c('breast_erpos','breast_erneg','endo','ovary','ovary_hgsoc','prostate')
for (pheno in phenos){
  pwas<-fread('/home/isabelle/FUSION/FUSION_PWAS_'%&%pheno%&%'_FDR.txt',stringsAsFactors=F)
  Phenotype<-rep_len(pheno,length.out=(nrow(pwas)))%>%as.data.frame()
  pwas$Phenotype<-Phenotype
  combopwas<-rbind(combopwas,pwas)
}
fwrite(combopwas,'/home/isabelle/FUSION/FUSION_PWAS_combined_FDR.txt',sep='\t',quote=F,row.names=F)
