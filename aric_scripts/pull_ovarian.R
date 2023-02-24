#for any fdr-sig protein, pull that protein's ovarian pvalue
library(dplyr)
library(data.table)

aric<-fread("~/FUSION/FUSION_PWAS_combined_FDR.txt",stringsAsFactors=F)
colnames(aric)
#[1] "PANEL"        "FILE"         "ID"           "CHR"          "P0"          
#[6] "P1"           "HSQ"          "BEST.GWAS.ID" "BEST.GWAS.Z"  "EQTL.ID"     
#[11] "EQTL.R2"      "EQTL.Z"       "EQTL.GWAS.Z"  "NSNP"         "NWGT"        
#[16] "MODEL"        "MODELCV.R2"   "MODELCV.PV"   "TWAS.Z"       "TWAS.P"      
#[21] "fdr_pvalue"   "Phenotype"   

fdrsig<-aric%>%filter(fdr_pvalue<0.05)
key<-fdrsig%>%pull(FILE)%>%unique()

ovary<-aric%>%filter(Phenotype=='ovary',FILE%in%key)
ovaryp<-ovary%>%select(FILE,TWAS.P,TWAS.Z)
colnames(ovaryp)<-c('FILE','ovary_pvalue','ovary_zscore')
fdrsig1<-left_join(fdrsig, ovaryp, by=c('FILE'))

fwrite(fdrsig1, "~/FUSION/ARIC_FDR0.05_with-ovary.txt",row.names=F,quote=F,sep="\t")

endo<-aric%>%filter(Phenotype=='endo',FILE%in%key)
endop<-endo%>%select(FILE,TWAS.P,TWAS.Z)
colnames(endop)<-c('FILE','endo_pvalue','endo_zscore')
fdrsig2<-left_join(fdrsig, endop, by=c('FILE'))

fwrite(fdrsig2, "~/FUSION/ARIC_FDR0.05_with-endo.txt",row.names=F,quote=F,sep="\t")
