#!/bin/bash

#Usage: ./subset_eRP_by_chr.sh gencode.vM4.genes.gtf
#Usage: ./subset_eRP_by_chr.sh mouse_allChr_wTAD_weRP_wAllComp_071819.txt

for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX
do
  awk -v var=$CHR '{if ($1 == var){print}}' $1 > ${CHR}_$1
done
