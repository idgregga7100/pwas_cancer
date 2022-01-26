#!/bin/bash
#pull lines with those genes from that table in those slides

for gene in KDELC2 GSTM1 RSPO3 PLG ARL3 MMP10 ABO
do
 grep $gene combined_rho_INTERVAL.csv > pulled_$gene.csv
done

cat pulled_KDELC2.csv pulled_GSTM1.csv pulled_RSPO3.csv pulled_PLG.csv pulled_ARL3.csv pulled_MMP10.csv pulled_ABO.csv > pulled_results_rho.csv
