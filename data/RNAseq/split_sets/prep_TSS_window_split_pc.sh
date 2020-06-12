#!/bin/bash

#prepare proportion of states in 75kbp window around train, test, and ref TSSs of just protein coding genes
python window_around_TSS.py rnaTPM_train_pc.txt train_TSS_window_pc.txt 60000
python window_around_TSS.py rnaTPM_test_pc.txt test_TSS_window_pc.txt 60000
python window_around_TSS.py rnaTPM_ref_pc.txt ref_TSS_window_pc.txt 60000

date; time bedtools intersect -a train_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > trainTSS_windows_int_all.bed
date; time bedtools intersect -a train_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > trainTSS_windows_v_all.bed

date; time bedtools intersect -a test_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > testTSS_windows_int_all.bed
date; time bedtools intersect -a test_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > testTSS_windows_v_all.bed

date; time bedtools intersect -a ref_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > refTSS_windows_int_all.bed
date; time bedtools intersect -a ref_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > refTSS_windows_v_all.bed

IDEAS=/project/vision/Data/IDEAS/ideasVisionV20p8Seg
BED=.bed
for CT in B Cd4 Cd8 Cfue Cfum Clp Cmp Er4 Eryad Eryfl G1e Gmp Hpc7 Imk Mk Lsk Mep Mon Neu Nk
do
  for type in ref train test
  do
    date; time awk -v var=$IDEAS$CT$BED '{if ($4 == var) {print}}' ${type}TSS_windows_int_all.bed > ${type}TSS_windows_int_$CT$BED
    date; time awk '{print $1, $2, $3, $6"_"$7"_"$8}' OFS='\t' ${type}TSS_windows_int_$CT$BED > ${type}TSS_windows_istate_$CT$BED
    date; time bedtools groupby -i ${type}TSS_windows_istate_$CT$BED -g 1,2,3 -c 4 -o collapse > ${type}TSS_windows_collapse_$CT$BED
  done
done
