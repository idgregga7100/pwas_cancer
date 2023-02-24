#doing fdr on meta and combining phenos
"%&%" = function(a,b) paste(a,b,sep="")

phenos<-c('breast','endo','ovary','prostate')
for (pheno in phenos){
  pwas<-fread('/home/isabelle/metal/METAL_ARIC_disc-rep_'%&%pheno%&%'1.tbl',stringsAsFactors = F)
  colnames(pwas)<-c('MarkerName','Allele1','Allele2','Weight','Zscore','Pvalue','Direction')
  fdr_pvalue<-p.adjust(pwas$Pvalue, method = c('fdr'))
  pwasfdr<-cbind(pwas,fdr_pvalue)
  fwrite(pwasfdr, '/home/isabelle/metal/METAL_ARIC_disc-rep_'%&%pheno%&%'_FDR.txt', sep = '\t', quote = F,row.names = F)
}

combopwas<-fread('/home/isabelle/metal/METAL_ARIC_disc-rep_breast_FDR.txt',stringsAsFactors=F)
combopwas<-mutate(combopwas, Phenotype='breast')

phenos<-c('endo','ovary','prostate')
for (pheno in phenos){
  pwas<-fread('/home/isabelle/metal/METAL_ARIC_disc-rep_'%&%pheno%&%'_FDR.txt',stringsAsFactors=F)
  pwas<-mutate(pwas, Phenotype=pheno)
  combopwas<-rbind(combopwas,pwas)
}
fwrite(combopwas,'/home/isabelle/metal/METAL_ARIC_disc-rep_combined_FDR.txt',sep='\t',quote=F,row.names=F)