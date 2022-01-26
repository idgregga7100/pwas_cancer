library(dplyr)
library(data.table)
"%&%" = function(a,b) paste(a,b,sep="")

models<-c("PCAIR_baseline_models_rho0.1_zpval0.05_sig_genes","PrediXcan_SomaID")

for (model in models){
	spred<-read.csv(model%&%"_all_pheno.csv")
	print(model)
	bonf<-5e-08
	#bonf <- 0.05/(2*7*nrow(spred))
	signif <- filter(spred,-log10(P) > -log10(bonf))
	print(bonf)
	signif <- unique(signif)  
	write.csv(signif,file=model%&%"_5e-8_all_pheno.csv",quote=F,row.names=F)
}

#spred<-read.csv('combined_rho_INTERVAL.csv')
#bonf<-5e-08
#bonf<- 0.05/(7*nrow(spred))
#signif<-filter(spred,-log10(P)>-log10(bonf))
#print(bonf)
#signif<-unique(signif)
#write.csv(signif,file="combined_rho_INTERVAL_5e-8.csv",quote=F,row.names=F)
