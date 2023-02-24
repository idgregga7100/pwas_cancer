#elyse's script /home/egeoffroy/topmed/SPrediXcan/scripts/Manhattan_all_pheno.R
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(qqman))
suppressPackageStartupMessages(library(wesanderson))
options(warn=-1)
"%&%" = function(a,b) paste(a,b,sep="")

#the BP chromosome file for INTERVAL was made from all the genes in the gencode file so in theory it should work for these too? copied into this directory/name changed
load_bp_chrom <- function(pop, type){
        print('/home/isabelle/FUSION/BP_Chrom_files/'%&%pop%&%'_'%&%type%&%'_BP_Chrome.csv')
        BP_Chrome <- read.table('/home/isabelle/FUSION/BP_Chrom_files/'%&%pop%&%'_'%&%type%&%'_BP_Chrome.csv', sep = ',', header = T)
        colnames(BP_Chrome) <- c('PROBE_ID', 'ID', 'chr', 'BP')

        BP_Chrome <- BP_Chrome %>%  
                transform(chr = str_replace(chr, "chr", "")) 
        BP_Chrome <- transform(BP_Chrome, chr=as.numeric(chr))
        BP_Chrome <- transform(BP_Chrome, BP=as.numeric(BP))
        colnames(BP_Chrome)<-c('PROBE_ID','ID','CHR','BP')
        #i know changing colnames twice is silly just let it happen
#        BP_Chrome$GENE <- gsub("\\..*","",BP_Chrome$PROBE_ID)
#will need to merge by gene_name/ID not ENSG label/PROBE_ID because these results don't have ENSG names
        return(BP_Chrome)
}

get_files <- function(model_type){
	files <- list.files('/home/isabelle/FUSION/fusion_twas-master/cancer_pwas', pattern = model_type, recursive = T, full.names=T)
	return(files)
}

format_files <- function(files, model_type){
  output<-data.frame()
	for(file in files){
		if(!str_detect(file, 'sig_genes') && !str_detect(file, 'png') && !str_detect(file, 'tiff')){
			pheno <- str_split(file, '/')[[1]][7]
			print(pheno)
			print(file)		
	  	if(str_detect(file, "ALL")){
        		model <- "ALL"
        		BP_Chrome <- load_bp_chrom('ALL', model_type)
      } 
	  	if(str_detect(file, 'CAU')){
	  		model <- 'CAU'
	  		BP_Chrome <- load_bp_chrom('CAU', model_type)
  		} 
  		if(str_detect(file, 'CHN')){
	  		model <- 'CHN'
	  		BP_Chrome <- load_bp_chrom('CHN', model_type)
	  	} 
	  	if(str_detect(file, 'HIS')){
	  		model <- 'HIS'
	  		BP_Chrome <- load_bp_chrom('HIS', model_type)
	  	}
	  	if(str_detect(file, 'AFA')){
	  		model <- 'AFA'
	  		BP_Chrome <- load_bp_chrom('AFA', model_type)
	  	}
	  	if(str_detect(file, 'FUSION')){
	  		model<- 'FUSION'
	  		BP_Chrome<-load_bp_chrom('FUSION', model_type)
	  	}

	  	S_Pred_file <- read.table(file, header = T,  sep = ',')      
#      		S_Pred_file$GENE <- gsub("\\..*","",S_Pred_file$gene) 
        #not gonna be able to merge with this so not needed i don't think
    	S_Pred_file <- S_Pred_file[-c( 9)]
    	names(S_Pred_file)[names(S_Pred_file) == 'TWAS.P'] <- 'P'

    	GWAS <- merge(S_Pred_file, BP_Chrome, by = c('ID','CHR'))      
#      		GWAS <- S_Pred_file
  		GWAS <- na.omit(GWAS)
  		names(GWAS)
#	  	colnames(GWAS)[16] <- "CHR"
      	
#      		GWAS <- GWAS %>%  #added by Jenny
#        		transform(CHR = str_replace(CHR, "chr", "")) 
#      		GWAS<- transform(GWAS, CHR=as.numeric(CHR)) 
	  	if(str_detect(pheno, 'unfiltered')){
	    	pheno <- str_split(pheno, 'unfiltered_')[[1]][2]
    		pheno <- str_replace(pheno, '.csv', '')
  		} 
	  	if(str_detect(pheno, '0.05')){
	    	pheno <- str_split(pheno, '0.05_')[[1]][2]
    		pheno <- str_replace(pheno, '.csv', '')
	  	}
	  	if(str_detect(pheno, 'PWAS')){
	    	pheno <- str_split(pheno, 'PWAS_')[[1]][2]
	    	pheno <- str_replace(pheno, '.csv', '')
	  	} 
      		
      GWAS$Phenotype <- rep(pheno, nrow(GWAS))
      GWAS$Model <- rep(model, nrow(GWAS))
      names(GWAS)

      output <- rbind(output, GWAS)
		}
	}		
  return(output)
}

make_manhattan <- function(output, model_type, num_pheno){
	axis.set <- output %>% 
  		group_by(CHR) %>% 
  		summarize(center = (max(BP) + min(BP)) / 2)
	ylim <- abs(floor(log10(min(output$P)))) + 2 
	if(ylim < 10){
		ylim <- 10
	}
	sig <- 5e-8

	# Prepare the dataset
	output <- output %>%  
 	 	# Compute chromosome size
 	 	group_by(CHR) %>% 
 	 	summarise(chr_len=max(BP)) %>% 
  
 	 	# Calculate cumulative position of each chromosome
 	 	mutate(tot=cumsum(chr_len)-chr_len) %>%
 	 	select(-chr_len) %>%
  
 	 	# Add this info to the initial dataset
  		left_join(output, ., by=c("CHR"="CHR")) %>%
  
  		# Add a cumulative position of each SNP
  		arrange(CHR, BP) %>%
  		mutate( BPcum=BP+tot) %>%
    
  		# Filter SNP to make the plot lighter
		filter(-log10(P)>0.5)

		# Prepare X axis
		axisdf <- output %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#		colorsss <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#276419", "#B15928", "#000000", "#8E0152")
#facet_wrap(~Model, nrow = 5)
		p <- ggplot(output, aes(x = BPcum, y = -log10(P), 
                          color = Phenotype, shape = Phenotype)) +
  			geom_point(alpha = 0.75) +
  			geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  			scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) + 
  			scale_shape_manual(values=rep(c(15,16,17,18,4,5,6,7,8,9,10,15,16,17), 2)) + scale_color_manual(values = wes_palette("Zissou1", 7, type = "continuous")) +
  			scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  			scale_size_continuous(range = c(0.5,3)) +
  			labs(title = model_type, x = NULL, 
       			y = "-log10(p)") + 
  			theme_minimal() +
  			theme( legend.position = 'bottom', legend.title = element_text(size = 12), legend.text = element_text(size = 12), strip.text = element_text(size=13),  
    			panel.border = element_blank(),
    			panel.grid.major.x = element_blank(),
    			panel.grid.minor.x = element_blank(),
    			axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  			) + ggsave(model_type%&%'_Manhattan_all_pheno_5e8.png',device='png', width = 12, height = 8)

}



	phenotypes <- c('breast','breast_erneg','breast_erpos','endo','ovary','ovary_hgsoc','prostate')
	pops <- c('FUSION')
	model_types <- c('PWAS')
#	output <- data.frame()

	for(model_type in model_types){
		print(model_type)
		files <- get_files(model_type)
		print(files)
#		output <- format_files(files, model_type, output)
		output <- format_files(files, model_type)
		print(unique(output$Phenotype))
		write.csv(output, '/home/isabelle/FUSION/FUSION'%&%model_type%&%'_all_pheno.csv', quote = F, row.names=F)
		
		make_manhattan(output, model_type, len(unique(output$Phenotype)))		
	}


