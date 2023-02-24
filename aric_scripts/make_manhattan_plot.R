#script to generate complex manhattan plot, taken from 01Manhattan2.R
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(qqman))
suppressPackageStartupMessages(library(wesanderson))
options(warn=-1)
"%&%" = function(a,b) paste(a,b,sep="")

#this creates an R function to generate the manplot, shouldn't need to edit any of this chunk
make_manhattan <- function(output, model_type, num_pheno){
  
  axis.set <- output %>% 
    group_by(CHR) %>% 
    summarize(center = (max(BP) + min(BP)) / 2)
  ylim <- abs(floor(log10(min(output$P)))) + 2 
  if(ylim < 10){
    ylim <- 10
  }
  sig <- 5.6e-6 #Dr. W changed
  
  # Prepare the dataset
  output <- output %>%  
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
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
    geom_point(alpha = 0.75, size=4) + #Dr. W added size
    geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) + 
    scale_shape_manual(values=rep(c(15,16,17,18,4,5,6,7,8,9,10,15,16,17), 2)) + scale_color_manual(values = wes_palette("Zissou1", 28, type = "continuous")) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(title = "Discovery PWAS", x = NULL, 
         y = "-log10(p)") + 
    theme_minimal() + #Dr. W changed font sizes below
    theme( legend.position = 'bottom', legend.title = element_text(size = 14), legend.text = element_text(size = 14), strip.text = element_text(size=14),  
           panel.border = element_blank(),
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
           axis.text.y = element_text(size = 12, vjust = 0.5)
    ) #+ ggsave(model_type%&%'_Manhattan_all_pheno_5e8.png',device='png', width = 12, height = 8)
  return(p) #Dr. W added this to return the plot from the function
}

phenotypes <- c('LDL_cholesterol', 'Fasting_glucose', 'Fasting_insulin', 'Coffee', 'BMI', 'C-reactive', 'Chronic_kidney', 'Diabetes', 'Diastolic_blood', 'Hemoglobin', 'Glomerular','Mean_corpuscular_hemoglobin', 'HDL_cholesterol', 'Total_cholesterol', 'Triglyceride','Systolic_blood', 'Smoking', 'WBC', 'Height','Waist-hip50', 'Waist-hip51', 'Waist-hip52', 'End_renal', 'QRS_duration', 'QT_interval', 'PR_interval', 'Hypertension', 'Platelet')
pops <- c('AFA', 'HIS', 'CAU', 'ALL', 'CHN')
model_type <- c('PCAIRbaseline_PAV_filtered_rho0.1_zpval0.05', 'PCAIRbaseline_PAV_filtered_unfiltered', 'PCAIRbaseline_rho0.1_zpval0.05', 'PCAIRbaseline_unfiltered')

#made the output with the format_SPrediXcan_output.R script (used to be one script but i split it in two)
output<-read.csv('/path/to/formatted/output/file')

#Dr. W added this and played with the dimensions and size parameters in the make_manhattan function
#to get the plot looking decent (just trial and error)
manplot <- make_manhattan(output, model_type, len(unique(output$Phenotype)))		
tiff("ARIC_PWAS_manplot.tiff", width = 20, height = 12, units = 'cm', res = 600, compression = 'lzw')
manplot
dev.off()