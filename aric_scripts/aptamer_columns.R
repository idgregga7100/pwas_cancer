#SeqId<-str_sub(condensed$FILE, start=66)
"%&%" = function(a,b) paste(a,b,sep="")

phenos<-c('breast','endo','ovary','prostate')
for (pheno in phenos){
#discovery data, ARIC models
data<-data<-fread('FUSION/FUSION_PWAS_'%&%pheno%&%'_FDR.txt')
SeqId<-str_sub(data$FILE, start=67, end=-10)%>%as.data.frame()
colnames(SeqId)<-c('SeqId')
newdata<-cbind(SeqId,data)
fwrite(newdata, 'FUSION/FUSION_PWAS_'%&%pheno%&%'_FDR_aptamers.txt',quote=F,row.names = F,sep='\t')
}

phenos<-c('breast','endom','ovary','prostate')
for (pheno in phenos){
#replication data, ARIC models
data<-data<-fread('pwas_cancer/replication/'%&%pheno%&%'_FUSION_PWAS.csv')
SeqId<-str_sub(data$FILE, start=66, end=-10)%>%as.data.frame()
colnames(SeqId)<-c('SeqId')
newdata<-cbind(SeqId,data)
fwrite(newdata, 'pwas_cancer/replication/'%&%pheno%&%'_FUSION_PWAS_aptamers.txt',quote=F,row.names = F,sep='\t')
}
