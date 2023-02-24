#combine the separate chromosome files into single phenotype files
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

phenos<-c('breast_erneg','breast_erpos','breast','endo','ovary','ovary_hgsoc','prostate')
#phenos<-c('breast_erneg')

for (pheno in phenos){
  file1<-fread('/home/isabelle/FUSION/fusion_twas-master/cancer_pwas/'%&%pheno%&%'_hg38_chr1.dat',header=T,stringsAsFactors=F)
  for (n in 2:22){
    file2<-fread('/home/isabelle/FUSION/fusion_twas-master/cancer_pwas/'%&%pheno%&%'_hg38_chr'%&%n%&%'.dat',header=T,stringsAsFactors=F)
    print(pheno%&%' chr'%&%n%&%' with '%&%nrow(file2)%&%' genes')
    file1<-rbind(file1,file2)
  }
  mhc<-fread('/home/isabelle/FUSION/fusion_twas-master/cancer_pwas/'%&%pheno%&%'_hg38_chr6.dat.MHC',header=T,stringsAsFactors=F)
  print(pheno%&%'MHC with '%&%nrow(mhc)%&%' genes')
  file1<-rbind(file1,mhc)
  print(pheno%&%' with '%&%nrow(file1)%&%' total genes')
  write.csv(file1,'/home/isabelle/FUSION/fusion_twas-master/cancer_pwas/'%&%pheno%&%'_FUSION_PWAS.csv',quote=F,row.names=F)
}