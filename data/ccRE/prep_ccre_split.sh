#!/bin/bash

#prepare proportion of states of train, test, and ref CREs
date; time bedtools intersect -a ../../usevision/train_cre_uniq.bed -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames >train_cre_int_all.bed
date; time bedtools intersect -a ../../usevision/train_cre_uniq.bed -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > train_cre_v_all.bed
wc -l train_cre_v_all.bed

date; time bedtools intersect -a ../../usevision/test_cre_uniq.bed -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > test_cre_int_all.bed
date; time bedtools intersect -a ../../usevision/test_cre_uniq.bed -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > test_cre_v_all.bed
wc -l test_cre_v_all.bed

date; time bedtools intersect -a ../../usevision/ref_cre_uniq.bed  -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > ref_cre_int_all.bed
date; time bedtools intersect -a ../../usevision/ref_cre_uniq.bed  -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > ref_cre_v_all.bed
wc -l ref_cre_v_all.bed

IDEAS=/project/vision/Data/IDEAS/ideasVisionV20p8Seg
BED=.bed
for CT in B Cd4 Cd8 Cfue Cfum Clp Cmp Er4 Eryad Eryfl G1e Gmp Hpc7 Imk Mk Lsk Mep Mon Neu Nk
do
  for type in ref train test
  do
    date; time awk -v var=$IDEAS$CT$BED '{if ($4 == var) {print}}' ${type}_cre_int_all.bed > ${type}_cre_int_$CT$BED
    date; time awk '{print $1, $2, $3, $6"_"$7"_"$8}' OFS='\t' ${type}_cre_int_$CT$BED > ${type}_cre_istate_$CT$BED
    date; time bedtools groupby -i ${type}_cre_istate_$CT$BED -g 1,2,3 -c 4 -o collapse > ${type}_cre_collapse_$CT$BED
  done
done

date; time python ccre_prop_in_state.py train_cre_collapse_{}.bed 179314 train_cre_state_prop.npz
date; time python ccre_prop_in_state.py test_cre_collapse_{}.bed 19293 test_cre_state_prop.npz
date; time python ccre_prop_in_state.py ref_cre_collapse_{}.bed 6412 ref_cre_state_prop.npz
