#!/bin/bash

for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX
do
  date; time ./find_protein_coding.py ${CHR}_mouse_allChr_wTAD_weRP_wAllComp_071819.txt ${CHR}_gencode.vM4.genes.gtf ${CHR}_pc_eRP_allInfo.txt
  cat ${CHR}_pc_eRP_allInfo.txt >> allCHR_pc_eRP_allInfo.txt
done
