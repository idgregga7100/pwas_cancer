#to run metal: two allele cols, which we decided can be like all A all G bcuz they don't matter for us
#also a common identifier column like isntead of a snp column: using ARIC seqids that i lined up in a huge key
#also also, i am making the decision that we're only running things that exist in ARIC bcuz those are the only results we're looking at anyway

#so i gotta add these cols
#files: discovery ARIC, INTERVAL, TOPMed; replication ARIC, TOPMed; within the four phenos (bcuz rep only has 4)
"%&%" = function(a,b) paste(a,b,sep="")
key<-fread('huge_key.tsv')

phenos<-c('breast','endo','ovary','prostate')
for (pheno in phenos){
aric<-fread('FUSION/FUSION_PWAS_'%&%pheno%&%'_FDR_aptamers.txt')
aricnew<-mutate(aric, A1='A')%>%mutate(A2='G')
fwrite(aricnew, 'metal/ARIC_discovery_'%&%pheno%&%'.txt',quote=F,row.names=F,sep='\t')

interval<-fread('topmed/proteome/Schlidkraut/INTERVAL_PrediXcan_'%&%pheno%&%'_FDR.txt')
intervalnew<-mutate(interval, A1='A')%>%mutate(A2='G')
intkey<-select(key, SeqId.x, gene)
intervalfinal<-left_join(intervalnew, intkey, by='gene')%>%filter(!is.na(SeqId.x))
fwrite(intervalfinal, 'metal/INTERVAL_discovery_'%&%pheno%&%'.txt',quote=F,row.names=F,sep='\t')

topmed<-fread('topmed/proteome/Schlidkraut/ALL_baseline_rhozpval_'%&%pheno%&%'_FDR.txt')
topmednew<-mutate(topmed, A1='A')%>%mutate(A2='G')
topkey<-select(key, SeqId.x, gene_name)
topmedfinal<-left_join(topmednew, topkey, by='gene_name')%>%filter(!is.na(SeqId.x))
fwrite(topmedfinal, 'metal/TOPMed_discovery_'%&%pheno%&%'.txt',quote=F,row.names=F,sep='\t')
}

phenos<-c('breast','endom','ovary','prostate')
for (pheno in phenos){
  aric<-fread('pwas_cancer/replication/'%&%pheno%&%'_FUSION_PWAS_aptamers.txt')
  aricnew<-mutate(aric, A1='A')%>%mutate(A2='G')
  fwrite(aricnew, 'metal/ARIC_replication_'%&%pheno%&%'.txt',quote=F,row.names=F,sep='\t')
  
  topmed<-fread('pwas_cancer/replication/ALL_baseline_rhozpval_'%&%pheno%&%'_FDR.txt')
  topmednew<-mutate(topmed, A1='A')%>%mutate(A2='G')
  topkey<-select(key, SeqId.x, gene_name)
  topmedfinal<-left_join(topmednew, topkey, by='gene_name')%>%filter(!is.na(SeqId.x))
  fwrite(topmedfinal, 'metal/TOPMed_replication_'%&%pheno%&%'.txt',quote=F,row.names=F,sep='\t')
}
