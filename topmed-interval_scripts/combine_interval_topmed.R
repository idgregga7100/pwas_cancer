#combine INTERVAL and TOPMed results so they're sire by side in one spreadsheet
#rho0.1 filtered models
library(stringr)
library(dplyr)
library(data.table)

topmedfile<-read.csv('/home/isabelle/topmed/proteome/Schlidkraut/PCAIR_baseline_models_PAV_filtered_rho0.1_zpval0.05_sig_genes_all_pheno.csv',stringsAsFactors=F)
#GENE,gene_name,gene,zscore,effect_size,P,var_g,pred_perf_r2,pred_perf_pval,n_snps_used,n_snps_in_cov,n_snps_in_model,PROBE_ID,CHR,BP,Phenotype,Model
#SL000003_ENSG00000214274,SL000003,SL000003_ENSG00000214274.9,-0.284266084537825,-0.00555382475477968,0.776206475087561,0.244622580794782,0.326411490895927,3.11772334679812e-91,52,52,52,SL000003_ENSG00000214274.9,14,20110251,breast_erneg,ALL
intervalfile<-read.csv('/home/isabelle/topmed/proteome/Schlidkraut/PrediXcan_sig_genes_all_pheno.csv',stringsAsFactors=F)
#gene_name,gene,zscore,effect_size,P,var_g,pred_perf_r2,pred_perf_pval,n_snps_used,n_snps_in_cov,n_snps_in_model,best_gwas_p,largest_weight,PROBE_ID,CHR,BP,Phenotype,Model
#ABO,ABO.9253.52.3,0.792206057471917,0.00719808267479617,0.428240537142174,0.707583025518869,0.745665597684641,0,56,61,61,0.0129299982233509,0.239830172187927,ENSG00000175164.16,9,133233278,breast_erneg,INTERVAL

topmed<-select(topmedfile, PROBE_ID)
id<-data.frame()[1:nrow(topmed),]
for(gene in topmed){
  gene<-as.character(gene)
  gene<-str_sub(gene,start=10,end=24)
  id<-cbind.data.frame(id,gene)
}
rownames(id)<-c(1:nrow(id))
colnames(id)<-c('ENSG_id')
topmedfile<-cbind(topmedfile,id)
#this adds a key of SL000000_ENSG00000000000 turned into ENSG00000000000 so that INTERVAL can be searched, cuts off .00

#ensg_labels<-pull(topmedfile, ENSG_id)
#ensg_labels<-as.character(ensg_labels)
#ensg_labels<-unique(ensg_labels)
#list of labels without any repeats except if doing full_join then don't actually need this
#print(ensg_labels)

topmedfileALL<-topmedfile%>%filter(Model=="ALL")%>%select(-gene,-Model)
colnames(topmedfileALL)<-c('GENE','gene_name','ALL_zscore','ALL_effect_size','ALL_P','ALL_var_g','ALL_pred_perf_r2','ALL_pred_perf_pval','ALL_n_snps_used','ALL_n_snps_in_cov','ALL_n_snps_in_model','PROBE_ID','CHR','BP','Phenotype','ENSG_id')
topmedfileCAU<-topmedfile%>%filter(Model=="CAU")%>%select(-gene,-Model)
colnames(topmedfileCAU)<-c('GENE','gene_name','CAU_zscore','CAU_effect_size','CAU_P','CAU_var_g','CAU_pred_perf_r2','CAU_pred_perf_pval','CAU_n_snps_used','CAU_n_snps_in_cov','CAU_n_snps_in_model','PROBE_ID','CHR','BP','Phenotype','ENSG_id)

ALLCAU<-full_join(topmedfileALL,topmedfileCAU)
print(head(ALLCAU))

interval<-select(intervalfile, gene,gene_name,Phenotype,PROBE_ID,zscore,effect_size,P,var_g,pred_perf_r2)
colnames(interval)<-c('gene','gene_name_id','Phenotype','PROBE_ID','INTERVAL_zscore','INTERVAL_effect_size','INTERVAL_P','INTERVAL_var_g','INTERVAL_pred_perf_r2')

genes<-select(interval, PROBE_ID)
id<-data.frame()[1:nrow(genes),]
for(gene in genes){
  gene<-as.character(gene)
  gene<-str_sub(gene,start=1,end=15)
  id<-cbind.data.frame(id,gene)
}
rownames(id)<-c(1:nrow(id))
colnames(id)<-c('ENSG_id')
interval<-cbind(interval,id)
#cuts off .00 from ENSG00000000000

combined<-full_join(ALLCAU,interval)
print(head(combined))
#write.csv(combined,'combo_INTERVAL_PAV_filtered_rho0.1_zpval0.05.csv',row.names=F,quote=F)
