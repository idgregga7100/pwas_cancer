#!/bin/bash

#Schlidkraut cancer colab data
#use Ryan's models, baseline
#/home/isabelle/topmed/proteome/Schlidkraut

p="CAU"
#p2="PAV_filtered_rho0.1_zpval0.05"
#p2="PAV_filtered_unfiltered"
#p2="rho0.1_zpval0.05"
p2="unfiltered"

FILES=/home/wheelerlab3/Data/Schlidkraut_BC_OC_collab/*_hg38.txt
for f in $FILES
do
  /usr/local/bin/MetaXcan_software/SPrediXcan.py \
  --model_db_path /home/ryan/TOPMed_Proteome/dbs_out/${p}_PCAIR_baseline_models_${p2}.db \
  --covariance /home/ryan/TOPMed_Proteome/dbs_out/${p}_PCAIR_baseline_models_${p2}_covariances.txt.gz \
  --gwas_file $f \
  --snp_column SNP_hg38 \
  --effect_allele_column A1 \
  --non_effect_allele_column A2 \
  --beta_column BETA \
  --se_column SE \
  --pvalue_column P \
  --output_file /home/isabelle/topmed/proteome/Schlidkraut/SPred_out/${p}_PCAIR_baseline_models_${p2}_${f:48:(-9)}.csv \
  --keep_non_rsid
done


#string manipulation in BASH https://www.tldp.org/LDP/abs/html/string-manipulation.html
#used to rm .txt.gz from output file name
#also the path from the file name
