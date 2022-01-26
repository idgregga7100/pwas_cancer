suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(qqman))
suppressPackageStartupMessages(library(wesanderson))
options(warn=-1)
#just trying to put all the predixcan interval phenotypes together with chunks from the big Manhattan_all_pheno.R script

get_files <- function(model_type){
        files <- list.files('/home/isabelle/topmed/proteome/Schlidkraut/INTERVAL_predixcan', pattern = model_type, recursive = T, full.names=T)
        return(files)
}

add_pheno_model <- function(files, output){
	for (file in files){
		pheno <- str_split(file, '/')[[1]][8]
		model<- 'INTERVAL'
		if(str_detect(pheno, 'PrediXcan')){
                	pheno <- str_split(pheno, 'PrediXcan_')[[1]][2]
                	pheno <- str_replace(pheno, '.csv', '')
               	}
        	GWAS <- read.table(file, header = T,  sep = ',')
		GWAS <- GWAS[-c( 9)]
                names(GWAS)[names(GWAS) == 'pvalue'] <- 'P'
		GWAS$Phenotype <- rep(pheno, nrow(GWAS))
        	GWAS$Model <- rep(model, nrow(GWAS))
		output<-rbind(output,GWAS)
	}
	return(output)
}

model_types <- c('PrediXcan')
output <- data.frame()

for(model_type in model_types){
        print(model_type)
        files <- get_files(model_type)
        print(files)
        output <- add_pheno_model(files, output)

	write.csv(output, paste('/home/isabelle/topmed/proteome/Schlidkraut/', model_type, '_SomaID_all_pheno.csv', sep = ''), quote = F, row.names=F)
}
