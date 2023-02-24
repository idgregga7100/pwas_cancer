#!/bin/bash

#run FUSION with sumstats and ARIC weights from Chatterjee lab, runs by chromosome

FILES=/home/isabelle/summary_stats/Rashkin_cancer/*_hg38_FUSION.txt

for f in ${FILES}
do

for n in {1..22}
do

Rscript /home/isabelle/FUSION/fusion_twas-master/FUSION.assoc_test.R \
--sumstats ${f} \
--weights /home/isabelle/Chatterjee_pwas/PWAS_EA/Plasma_Protein_EA_hg38.pos \
--weights_dir /home/isabelle/Chatterjee_pwas/PWAS_EA/Plasma_Protein_weights_EA \
--ref_ld_chr /home/isabelle/Chatterjee_pwas/LDref/EUR/chr \
--chr ${n} \
--out /home/isabelle/pwas_cancer/replication/${f:44:(-16)}_chr${n}.dat

done

done

#string manipulation to remove path etc/isolate pheno from file name
