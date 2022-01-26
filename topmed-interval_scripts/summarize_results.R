library(dplyr)
library(data.table)
library(stringr)

#replace CAU with EUR
file<-read.csv('pulled_results_combined_rho_INTERVAL.csv',stringsAsFactors=F) 
#combined_rho_INTERVAL.csv combined_PAV_rho_INTERVAL.csv pulled_results_combined_PAV_rho_INTERVAL.csv pulled_results_combined_rho_INTERVAL.csv
filecols<-colnames(file)
filecolsEUR<-filecols%>%str_replace("CAU","EUR")
colnames(file)<-c(filecolsEUR)

#get bonferroni significant for topmed--2 pops, 7 phenos, once again am i doing this right idk
#bonf <- 0.05/(7*2*nrow(file))
bonf <- 0.05
signifor <- filter(file,-log10(ALL_P) > -log10(bonf) | -log10(EUR_P) > -log10(bonf))
signifand <- filter(file,-log10(ALL_P) > -log10(bonf) & -log10(EUR_P) > -log10(bonf))
signif<-rbind(signifand,signifor)%>%unique()
#idk if there's an and/or symbol but this will work right?
print(bonf)
write.csv(signif,'pulled_results_combined_rho_INTERVAL_p0.05.csv',row.names=F,quote=F)

