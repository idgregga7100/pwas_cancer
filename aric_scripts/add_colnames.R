#add the header whoops idk how to do it in bash
library(data.table)
file1<-read.csv('FUSIONPWAS_all_pheno.csv',stringsAsFactors=F)
cols<-colnames(file1)

file2<-read.csv('FUSION_PWAS_pulled_results.csv',stringsAsFactors=F)
colnames(file2)<-c(cols)

write.csv(file2,'FUSION_PWAS_pulled_results.csv',quote=F,row.names=F)