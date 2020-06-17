#!/bin/bash

#prepare proportion of states in 60kbp window around train, test, and ref TSSs of just protein coding genes
#the output will be chr window_start window_stop gene_name
#python window_around_TSS.py rnaTPM_train_pc.txt train_TSS_window_pc.txt 60000
#python window_around_TSS.py rnaTPM_test_pc.txt test_TSS_window_pc.txt 60000
#python window_around_TSS.py rnaTPM_ref_pc.txt ref_TSS_window_pc.txt 60000

#wc -l *_TSS_window_pc.txt

#date; time bedtools intersect -a train_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > trainTSS_windows_int_all.bed
#date; time bedtools intersect -a train_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > trainTSS_windows_v_all.bed

#date; time bedtools intersect -a test_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > testTSS_windows_int_all.bed
#date; time bedtools intersect -a test_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > testTSS_windows_v_all.bed

#date; time bedtools intersect -a ref_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > refTSS_windows_int_all.bed
#date; time bedtools intersect -a ref_TSS_window_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > refTSS_windows_v_all.bed

#wc -l *TSS_windows_int_all.bed
#wc -l *TSS_windows_v_all.bed

IDEAS=/project/vision/Data/IDEAS/ideasVisionV20p8Seg
BED=.bed
for CT in B Cd4 Cd8 Cfue Cfum Clp Cmp Er4 Eryad Eryfl G1e Gmp Hpc7 Imk Mk Lsk Mep Mon Neu Nk
do
  for type in ref train test
  do
    #awk -v var=$IDEAS$CT$BED '{if ($5 == var) {print}}' ${type}TSS_windows_int_all.bed > ${type}TSS_windows_int_$CT$BED
    #awk '{print $1, $2, $3, $4, $6"_"$7"_"$8}' OFS='\t' ${type}TSS_windows_int_$CT$BED > ${type}TSS_windows_istate_$CT$BED
    #bedtools groupby -i ${type}TSS_windows_istate_$CT$BED -g 1,2,3,4 -c 5 -o collapse > ${type}TSS_windows_collapse_$CT$BED
    #echo $type $CT
    #wc -l ${type}TSS_windows_*_$CT$BED
    awk '{print $1, $2, $3, $5}' OFS='\t' ${type}TSS_windows_collapse_$CT$BED > ${type}TSS_windows_fcollapse_$CT$BED #remove the gene name now so that the field numbers match in the downstream script
  done
done

#wc -l trainTSS_windows_collapse_*.bed #19003?
#wc -l testTSS_windows_collapse_*.bed #1966?
#wc -l refTSS_windows_collapse_*.bed #857?

date; time python ccre_prop_in_state.py trainTSS_windows_fcollapse_{}.bed 19012 trainTSS_window_state_prop_pc.npz
date; time python ccre_prop_in_state.py testTSS_windows_fcollapse_{}.bed 1967 testTSS_window_state_prop_pc.npz
date; time python ccre_prop_in_state.py refTSS_windows_fcollapse_{}.bed 857 refTSS_window_state_prop_pc.npz
