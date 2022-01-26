#add the header whoops idk how to do it in bash
library(data.table)
file1<-read.csv('combined_rho_INTERVAL.csv',stringsAsFactors=F)
cols<-colnames(file1)

file2<-read.csv('pulled_results_PAV_rho.csv',stringsAsFactors=F)
colnames(file2)<-c(cols)

write.csv(file2,'pulled_results_combined_PAV_rho_INTERVAL.csv',quote=F,row.names=F)

file3<-read.csv('pulled_results_rho.csv',stringsAsFactors=F)
colnames(file3)<-c(cols)

write.csv(file3,'pulled_results_combined_rho_INTERVAL.csv',quote=F,row.names=F)
