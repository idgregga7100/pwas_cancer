#need to add z score, fix colnames: SNP, A1, A2, Z
#trying rsid as SNP, snp column just named hg38?
library(stringr)
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

phenos<-c('breast_erneg','breast_erpos','breast','endo','ovary','ovary_hgsoc','prostate')

for (pheno in phenos){
#	gwas<-read.table('/home/wheelerlab3/Data/Schlidkraut_BC_OC_collab/'%&%pheno%&%'_hg38.txt',stringsAsFactors=FALSE,header=TRUE)
	gwas<-read.table('/home/isabelle/Chatterjee_pwas/Schlidkraut/sumstats/'%&%pheno%&%'_hg38_rsid.txt',stringsAsFactors=FALSE,header=TRUE)
#	gwas['Z']<-NA
#	gwas$Z <- gwas$BETA/gwas$SE
#	colnames(gwas)<-c('CHR','BP','SNP','hg38','A1','A2','FREQ','BETA','SE','P','Z')
#	gwas$Z<-str_replace(gwas$Z, 'NaN', 'NA')
#	gwas<-na.omit(gwas)
	gwas<-gwas%>%select(SNP,A1,A2,Z)
	print(head(gwas))
	fwrite(gwas, file='/home/isabelle/Chatterjee_pwas/Schlidkraut/sumstats/'%&%pheno%&%'_hg38_cols.txt', quote=FALSE,row.names=FALSE,sep=' ',append=FALSE,col.names=TRUE)
}
