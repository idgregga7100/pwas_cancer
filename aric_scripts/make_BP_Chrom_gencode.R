#read gencode .gtf file, extract attributes, make BP Chrom file for INTERVAL PrediXcan results
library(data.table)
library(dplyr)
genes<- fread("/home/isabelle/topmed/proteome/Schlidkraut/gencode.v38.annotation.gtf")
colnames(genes)<-c('seqname','source','feature','start','end','score','strand','frame','attributes')
genes<-genes%>%filter(feature=='gene')

#function to deal with the attributes field
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}

genes$gene_id <- unlist(lapply(genes$attributes, extract_attributes, "gene_id"))
genes$gene_name <- unlist(lapply(genes$attributes, extract_attributes, "gene_name"))

#BP Chrom files are gene gene_name chr BP
#using start as the BP because idk what else i'd do
bp_chrom<-genes%>%select(gene_id, gene_name, seqname, start)
colnames(bp_chrom)<-c("gene","gene_name","chr","BP")
write.csv(bp_chrom, "/home/isabelle/topmed/proteome/Schlidkraut/BP_Chrom_files/INTERVAL_PrediXcan_BP_Chrom.csv",quote=F,row.names=F)
#INTERVAL Pred output: will need to merge these by gene_name column, not by gene column
