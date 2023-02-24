#i don't trust my previous effort to get the overlapping proteins between the models so we're trying again for a better/simpler result
library(dplyr)
library(data.table)
library(stringr)

#i have this huge combined file this time! hell yeah!
aric_fusion<-fread('~/FUSION/Cancer_PrediXcan_with_FDR.txt',stringsAsFactors = FALSE)

goi<-aric_fusion%>%filter(!is.na(CHR))%>%select(FILE, INT_gene, ID)%>%filter(!is.na(INT_gene))%>%unique()

seqid<-str_split_fixed(goi$FILE, '/', n=7)
seqid<-seqid[,7]
SeqId<-str_sub(seqid, start=2,end=-10)
goi$FILE<-SeqId
colnames(goi)<-c('ARIC_name','INT_name','ID')

fwrite(goi, 'ARIC_tested_in_INTERVAL.txt',quote=F,row.names=F,sep='\t')

goi<-aric_fusion%>%filter(!is.na(CHR))%>%select(FILE, Aptamer, ALL_P, ID)%>%filter(!is.na(ALL_P))%>%select(-ALL_P)%>%unique()

seqid<-str_split_fixed(goi$FILE, '/', n=7)
seqid<-seqid[,7]
SeqId<-str_sub(seqid, start=2,end=-10)
goi$FILE<-SeqId
colnames(goi)<-c('ARIC_name','MESA_name','ID')

fwrite(goi, 'ARIC_tested_in_TOPMedMESA-ALL.txt',quote=F,row.names=F,sep='\t')

goi<-aric_fusion%>%filter(!is.na(CHR))%>%select(FILE, INT_gene, Aptamer, ALL_P, ID)
goi<-filter(goi, !is.na(INT_gene) & !is.na(ALL_P))%>%select(-ALL_P)%>%unique()

seqid<-str_split_fixed(goi$FILE, '/', n=7)
seqid<-seqid[,7]
SeqId<-str_sub(seqid, start=2,end=-10)
goi$FILE<-SeqId
colnames(goi)<-c('ARIC_name','INT_name','MESA_name','ID')

fwrite(goi, 'ARIC_tested_in_MESA-INT.txt',quote=F,row.names=F,sep='\t')
