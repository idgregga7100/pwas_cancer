#ggpairs(fusionfdr, columns='TWAS.Z', aes(color=Phenotype))+ggsave('FUSION_pheno_zscore_correlation.png',device='png')
library(data.table)
library(dplyr)
library(GGally)
"%&%" = function(a,b) paste(a,b,sep="")

phenos_add<-c('breast_erpos','breast','endo','ovary','ovary_hgsoc','prostate')

file1<-fread('/home/isabelle/FUSION/FUSION_PWAS_breast_erneg_FDR.txt',header=T,stringsAsFactors=F)
id_z_main<-select(file1, FILE, TWAS.Z, fdr_pvalue)
colnames(id_z_main)<-c('ID','Z_breast_erneg','fdr_p_breast_erneg')

for (pheno in phenos_add){
  file2<-fread('/home/isabelle/FUSION/FUSION_PWAS_'%&%pheno%&%'_FDR.txt',header=T,stringsAsFactors=F)
  id_z<-select(file2, FILE, TWAS.Z, fdr_pvalue)
  colnames(id_z)<-c('ID','Z_'%&%pheno,'fdr_p_'%&%pheno)
  id_z_main<-left_join(id_z_main, id_z, by='ID')
}

#print(head(id_z_main))

#unfiltered plot, made before including fdr columns
#ggpairs(id_z_main, columns=2:8)
#ggsave('FUSION_pheno_zscore_correlation2.png',device='png')

x<-c(0.05,0.1,0.2,0.3,0.4)
#n_proteins_all<-data.table(matrix(nrow=1,ncol=2))
#colnames(n_proteins_all)<-c('threshold','N')
for(y in x){
  filtered<-id_z_main%>%filter(fdr_p_breast_erneg < y | fdr_p_breast_erpos < y | fdr_p_breast < y | fdr_p_endo < y | fdr_p_ovary < y | fdr_p_ovary_hgsoc < y | fdr_p_prostate < y)
#  n_proteins<-data.table()
#  n_proteins$threshold<-y
#  n_proteins$N<-nrow(filtered)
#  print(n_proteins)
  #n_proteins_all<-rbind(n_proteins_all,n_proteins)
  just_Z<-select(filtered, contains('Z'))
  ggpairs(just_Z, lower=list(continuous=wrap('smooth',size=0.1)))
  ggsave('FUSION_zscore_correlation_any'%&%y%&%'_fitline_small.png',device='png',width=7, height=7)
}
#the loop won't work idfk

filtered<-id_z_main%>%filter(fdr_p_breast_erneg < 0.05 | fdr_p_breast_erpos < 0.05 | fdr_p_breast < 0.05 | fdr_p_endo < 0.05 | fdr_p_ovary < 0.05 | fdr_p_ovary_hgsoc < 0.05 | fdr_p_prostate < 0.05)
just_Z<-select(filtered, contains('Z'))
ggpairs(just_Z, lower=list(continuous=wrap('smooth',size=0.2)))
ggsave('FUSION_zscore_correlation_any0.05_fitline_small.png',device='png',width=7, height=7)


#fdrphenos<-c('fdr_p_breast_erneg','fdr_p_breast_erpos','fdr_p_breast','fdr_p_endo','fdr_p_ovary','fdr_p_ovary_hgsoc','fdr_p_prostate')
#for (indiv in fdrphenos){
#  print(indiv)
#  filtered<-id_z_main%>%filter(indiv < 0.1)
#  filtered_Z<-select(filtered, contains('Z'))
#  ggpairs(filtered_Z)
#  ggsave('FUSION_zscore_correlation_'%&%indiv%&%'0.1.png',device='png')
#}

