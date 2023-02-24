#oh my god getting all these aptamer names lined up and merged jfc gonna be the death of me

aric<-fread('FUSION/FUSION_PWAS_breast_FDR_aptamers.txt')
aricids<-select(aric, SeqId, ID)
digits<-str_sub(aricids$SeqId, start=7)%>%as.data.frame()
colnames(digits)<-c('digits')
aricids<-cbind(aricids,digits)

interval<-fread('topmed/proteome/Schlidkraut/INTERVAL_PrediXcan_breast_FDR.txt')
intervalids<-select(interval, gene, gene_name)
digits<-str_split_fixed(intervalids$gene, pattern='\\.',n=2)%>%as.data.frame()%>%select(V2)
digits<-str_replace_all(digits$V2, '\\.','_')%>%as.data.frame()
colnames(digits)<-c('digits')
intervalids<-cbind(intervalids,digits)

mesa<-fread('topmed/proteome/Schlidkraut/ALL_baseline_rhozpval_breast_FDR.txt')
mesaids<-select(mesa,gene_name)
key<-select(topmedaptamers, SeqId, SomaId, EntrezGeneSymbol)
mesaids<-left_join(mesaids,key, by=c('gene_name'='SomaId'))
digits<-str_replace_all(mesaids$SeqId, '-','_')%>%as.data.frame()
colnames(digits)<-c('digits')
mesaids<-cbind(mesaids,digits)

#totalkey<-merge(aricids,mesaids,all=T) #for apparently cursed evil reasons going by digits is NOT gonna work

aric2<-select(aricids, -digits)
interval2<-select(intervalids, -digits)
mesa2<-select(mesaids, -digits)
totalkey<-full_join(aric2,interval2, by=c('ID'='gene_name'))
totalkey<-full_join(totalkey,mesa2, by=c('ID'='EntrezGeneSymbol'))
#and this isn't perfect in any way shape or form so for my merging and joining sanity i'm just gonna fix the repeats and multiple aptamers by hand in docs