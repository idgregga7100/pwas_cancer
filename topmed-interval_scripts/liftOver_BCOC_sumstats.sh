#!/bin/bash
#retrieve hg19 to hg38 liftOver chain file
# find LiftOver files link at http://hgdownload.soe.ucsc.edu/downloads.html#human
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# make liftOver bed format input, awk pulls coordinates, tail removes header (first line)
awk '{print "chr"$1,$2,$2+1}' < breast.txt | tail -n+2 > cancer.bed

#run liftOver
liftOver cancer.bed hg19ToHg38.over.chain.gz cancer_hg38.bed cancer_unmapped.bed
#note, only need to do once because all sumstats files have the same SNPs

#I discovered endo.txt was in a different order, so I reordered the columns first to match the others
# don't run these commands again
#awk '{print $1, $2, $3, $5, $6, $4, $7, $8, $9}' endo.txt > o
#mv o endo.txt

#run python script to update coordinates and format sumstats for S-PrediXcan
#this is a bash for loop, I always google how to do it
for gwas in breast_erneg breast_erpos breast endo ovary_hgsoc ovary prostate
do
    python3 sumstats_hg38.py -s ${gwas}.txt -m cancer_hg38.bed -u cancer_unmapped.bed -o ${gwas}_hg38.txt
done
