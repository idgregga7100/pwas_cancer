#okay trying this again
#sticking topmed ALL and CAU results side by side, then adding INTERVAL

library(data.table)
library(dplyr)
library(stringr)

topmedfile<-read.csv('/home/isabelle/topmed/proteome/Schlidkraut/PCAIR_baseline_models_rho0.1_zpval0.05_sig_genes_all_pheno.csv',stringsAsFactors = FALSE)
#PCAIR_baseline_models_PAV_filtered_rho0.1_zpval0.05_sig_genes_all_pheno.csv
#PCAIR_baseline_models_rho0.1_zpval0.05_sig_genes_all_pheno.csv

ensg<-select(topmedfile, PROBE_ID)
id<-data.frame()[1:nrow(ensg),]
for(gene in ensg){
  gene<-as.character(gene)
  gene<-str_sub(gene,start=10,end=24)
  id<-cbind.data.frame(id,gene)
}
rownames(id)<-c(1:nrow(id))
colnames(id)<-c('ENSG_id')
topmedfile<-cbind(topmedfile,id)

topmedALL<-topmedfile%>%filter(Model=="ALL")%>%select(-gene,-Model)
colnames(topmedALL)<-c('GENE','gene_name','ALL_zscore','ALL_effect_size','ALL_P','ALL_var_g','ALL_pred_perf_r2','ALL_pred_perf_pval','ALL_n_snps_used','ALL_n_snps_in_cov','ALL_n_snps_in_model','PROBE_ID','CHR','BP','Phenotype','ENSG_id')
topmedCAU<-topmedfile%>%filter(Model=="CAU")%>%select(-gene,-Model)
colnames(topmedCAU)<-c('GENE','gene_name','CAU_zscore','CAU_effect_size','CAU_P','CAU_var_g','CAU_pred_perf_r2','CAU_pred_perf_pval','CAU_n_snps_used','CAU_n_snps_in_cov','CAU_n_snps_in_model','PROBE_ID','CHR','BP','Phenotype','ENSG_id')

topmed<-full_join(topmedALL,topmedCAU)
colnames(topmed)

intervalfile<-read.csv('/home/isabelle/topmed/proteome/Schlidkraut/PrediXcan_sig_genes_all_pheno.csv',stringsAsFactors = FALSE)

ensg<-select(intervalfile, PROBE_ID)
id<-data.frame()[1:nrow(ensg),]
for(gene in ensg){
  gene<-as.character(gene)
  gene<-str_sub(gene,start=1,end=15)
  id<-cbind.data.frame(id,gene)
}
rownames(id)<-c(1:nrow(id))
colnames(id)<-c('ENSG_id')
intervalfile<-cbind(intervalfile,id)

interval<-intervalfile%>%select(gene,Phenotype,zscore,effect_size,P,var_g,pred_perf_r2,ENSG_id)
colnames(interval)<-c('gene','Phenotype','INTERVAL_zscore','INTERVAL_effect_size','INTERVAL_P','INTERVAL_var_g','INTERVAL_pred_perf_r2','ENSG_id')

combo<-full_join(topmed,interval)
colnames(combo)

col_order <- c("GENE", "gene_name", "gene",'ENSG_id','Phenotype','ALL_zscore','ALL_effect_size','ALL_P','ALL_var_g','ALL_pred_perf_r2','ALL_pred_perf_pval','ALL_n_snps_used','ALL_n_snps_in_cov','ALL_n_snps_in_model','CAU_zscore','CAU_effect_size','CAU_P','CAU_var_g','CAU_pred_perf_r2','CAU_pred_perf_pval','CAU_n_snps_used','CAU_n_snps_in_cov','CAU_n_snps_in_model','INTERVAL_zscore','INTERVAL_effect_size','INTERVAL_P','INTERVAL_var_g','INTERVAL_pred_perf_r2','PROBE_ID','CHR','BP')
combo_ordered <- combo[, col_order]
colnames(combo_ordered)

write.csv(combo_ordered,'combined_rho_INTERVAL.csv',row.names=F,quote=F)
