#DO NOT RUN THIS AS-IS FROM TERMINAL COMMAND LINE gotta fidget with it, she's messy

library(dplyr)
library(stringr)
library(data.table)

interval<-read.csv('~/topmed/proteome/Schlidkraut/PrediXcan_sig_genes_all_pheno.csv',stringsAsFactors = FALSE) #PrediXcan_sig_genes_all_pheno.csv
#topmed_pavrho<-read.csv('~/topmed/proteome/Schlidkraut/PCAIR_baseline_models_PAV_filtered_rho0.1_zpval0.05_bonferroni_all_pheno.csv',stringsAsFactors = FALSE) #PCAIR_baseline_models_PAV_filtered_rho0.1_zpval0.05_sig_genes_all_pheno.csv 
topmed_rho<-read.csv('~/topmed/proteome/Schlidkraut/PCAIR_baseline_models_rho0.1_zpval0.05_sig_genes_all_pheno.csv',stringsAsFactors = FALSE) #PCAIR_baseline_models_rho0.1_zpval0.05_sig_genes_all_pheno.csv 
aric_fusion<-read.csv('~/FUSION/FUSION_PWAS_5.6e-06_all_pheno.csv',stringsAsFactors = FALSE)

#aric_fusion, interval
#ENSG00000134184.13
long_ensg<-pull(aric_fusion,PROBE_ID)
long_ensg<-as.character(long_ensg)
ensg<-long_ensg%>%str_sub(start=1,end=15)
ensg<-as.data.frame(ensg)
aric_fusion<-cbind(ensg,aric_fusion)


#topmed_pavrho, topmed_rho
#SL003739_ENSG00000243509.6
sl_ensg<-pull(topmed_rho,PROBE_ID)
sl_ensg<-as.character(sl_ensg)
ensg<-sl_ensg%>%str_sub(start=10,end=24)
ensg<-as.data.frame(ensg)
topmed_rho<-cbind(ensg,topmed_rho)


aric_fusion_genes_list<-pull(aric_fusion,ensg)
aric_fusion_genes_list<-as.character(aric_fusion_genes_list)
aric_fusion_genes_list<-unique(aric_fusion_genes_list)

interval_overlap<-filter(interval,ensg%in%aric_fusion_genes_list)
#topmed_pavrho_overlap<-filter(topmed_pavrho,ensg%in%aric_fusion_genes_list)
topmed_rho_overlap<-filter(topmed_rho,ensg%in%aric_fusion_genes_list)

write.csv(aric_fusion,'FUSION_PWAS_5.6e-06_all_pheno.csv',row.names = F,quote = F)
#now CAN CLEAR ALL THE STUFF

library(data.table)
library(dplyr)
library(tidyr)

fusion<-read.csv('FUSION_PWAS_5.6e-06_all_pheno.csv',stringsAsFactors = F)
fusion_interval<-read.csv('FUSION-INTERVAL_5.6e-06_overlap.csv',stringsAsFactors = F)
#fusion_topmedpavrho<-read.csv('FUSION-TOPMedPAVrho_bonf_overlap.csv',stringsAsFactors = F)
fusion_topmedrho<-read.csv('FUSION-TOPMed_5.6e-06_overlap.csv',stringsAsFactors = F)

#top_overlap<-inner_join(fusion_topmedpavrho,fusion_topmedrho,by=c('ensg','gene_name','Model'))
#toprho_interval_overlap<-inner_join(fusion_topmedrho,fusion_interval,by='ensg')

#int colnames
int_colnames<-colnames(fusion_interval)
int_string<-rep_len(c('INT_'),length(int_colnames))
int_colnames<-str_c(int_string,int_colnames)
colnames(fusion_interval)<-int_colnames

#need to split ALL/CAU, then colnames
ALL<-filter(fusion_topmedrho,Model=='ALL')
EUR<-filter(fusion_topmedrho,Model=='CAU')
ALL<-select(ALL,-Model)
EUR<-select(EUR,-Model)
ALL_colnames<-colnames(ALL)
EUR_colnames<-colnames(EUR)
ALL_string<-rep_len(c('ALL_'),length(ALL_colnames))
EUR_string<-rep_len(c('EUR_'),length(EUR_colnames))
ALL_colnames<-str_c(ALL_string,ALL_colnames)
colnames(ALL)<-ALL_colnames
EUR_colnames<-str_c(EUR_string,EUR_colnames)
colnames(EUR)<-EUR_colnames

topmed<-full_join(ALL,EUR, by= c('ALL_ensg'='EUR_ensg', 'ALL_gene_name'='EUR_gene_name', 'ALL_Phenotype'='EUR_Phenotype'))

replicates<-full_join(fusion_interval,topmed, by=c('INT_ensg'='ALL_ensg','INT_Phenotype'='ALL_Phenotype'))

aric_replicated<-left_join(fusion,replicates, by=c('ensg'='INT_ensg', 'ID'='INT_gene_name','Phenotype'='INT_Phenotype'))

write.csv(aric_replicated,'ARIC_5.6e-06_replicated.csv',row.names = F,quote = F)                      
#somehow this is still wrong idk whatever