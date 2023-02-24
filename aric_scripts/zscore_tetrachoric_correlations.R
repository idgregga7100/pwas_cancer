"%&%" = function(a,b) paste(a,b,sep="")
library(dplyr)
library(data.table)
library(psych)
library(stringr)

#create data.table of aric fdr and z for each pheno (1318 proteins 7 phenos)
phenos_add<-c('breast_erpos','breast','endo','ovary','ovary_hgsoc','prostate')

file1<-fread('/home/isabelle/topmed/proteome/Schlidkraut/CAU_baseline_rhozpval_breast_erneg_FDR.txt',header=T,stringsAsFactors=F)
id_z_main<-select(file1, gene, zscore, fdr_pvalue)
colnames(id_z_main)<-c('ID','Z_breast_erneg','fdr_p_breast_erneg')

for (pheno in phenos_add){
  file2<-fread('/home/isabelle/topmed/proteome/Schlidkraut/CAU_baseline_rhozpval_'%&%pheno%&%'_FDR.txt',header=T,stringsAsFactors=F)
  id_z<-select(file2, gene, zscore, fdr_pvalue)
  colnames(id_z)<-c('ID','Z_'%&%pheno,'fdr_p_'%&%pheno)
  id_z_main<-left_join(id_z_main, id_z, by='ID')
}

#i was gonna try to use stringr to extract seqids from the awful evil ID column but for whatever reason it hates me so i'm not gonna do that and i don't think it matters for this anyway ugh
just_Zall<-select(id_z_main, contains('Z'))
tetramatrix<-data.frame(matrix(nrow=nrow(just_Zall)))
for(z in just_Zall){
  tetra<-ifelse(z>0, 1, 0)%>%as.data.frame()
  tetramatrix<-cbind(tetramatrix,tetra)
}
tetramatrix<-tetramatrix[,2:8]
colnames(tetramatrix)<-c('breast_erneg','breast_erpos','breast','endo','ovary','ovary_hgsoc','prostate')
print(head(tetramatrix))
fwrite(tetramatrix, '/home/isabelle/FUSION/Zscoresunfiltered_for_tetra_TOPMedEUR.txt')

thresholds<-c(0.05,0.1,0.2,0.3,0.4)
for(n in thresholds){
  filtered<-id_z_main%>%filter(fdr_p_breast_erneg < n | fdr_p_breast_erpos < n | fdr_p_breast < n | fdr_p_endo < n | fdr_p_ovary < n | fdr_p_ovary_hgsoc < n | fdr_p_prostate < n)
  Z_filtered<-select(filtered, contains('Z'))
  print(nrow(Z_filtered))
  
  tetramatrix<-data.frame(matrix(nrow=nrow(Z_filtered)))
  for(z in Z_filtered){
    tetra<-ifelse(z>0, 1, 0)%>%as.data.frame()
    tetramatrix<-cbind(tetramatrix,tetra)
  }
  tetramatrix<-tetramatrix[,2:8]
  colnames(tetramatrix)<-c(n%&%'breast_erneg',n%&%'breast_erpos',n%&%'breast',n%&%'endo',n%&%'ovary',n%&%'ovary_hgsoc',n%&%'prostate')
  print(head(tetramatrix))
  fwrite(tetramatrix, '/home/isabelle/FUSION/Zscores'%&%n%&%'_for_tetra_TOPMedEUR.txt')
}

for(n in thresholds){
  zscores<-fread('/home/isabelle/FUSION/Zscores'%&%n%&%'_for_tetra_TOPMedEUR.txt', stringsAsFactors = F)
  print(n)
  output<-tetrachoric(zscores)
  table<-output$rho%>%as.data.table()
  rownames(table)<-c('breast_erneg','breast_erpos','breast','endo','ovary','ovary_hgsoc','prostate')
  fwrite(table, '/home/isabelle/FUSION/tetrachoric_'%&%n%&%'_TOPMedEUR.txt',row.names = T, quote = F)
}
zscores<-fread('/home/isabelle/FUSION/Zscoresunfiltered_for_tetra_TOPMedEUR.txt', stringsAsFactors = F)
output<-tetrachoric(zscores)
table<-output$rho%>%as.data.table()
rownames(table)<-c('breast_erneg','breast_erpos','breast','endo','ovary','ovary_hgsoc','prostate')
#output$rho = correlation table, the rownames are phenos but to make them actually show up gotta reset them
fwrite(table, '/home/isabelle/FUSION/tetrachoric_unfiltered_TOPMedEUR.txt',row.names = T, quote = F)
