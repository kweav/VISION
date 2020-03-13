#!/bin/bash

#prepare TPM files from count data
awk '{if ($1 != "chrM") {print}}' rnaHtseqCountsall_withcoordinates.0.txt > rnaHtseqCount_noM_withcoords.txt

echo chr$'\t'start$'\t'end$'\t'gene$'\t'gene_type$'\t'strand$'\t'Cfue$'\t'Cfum$'\t'Cmp$'\t'Eryfl$'\t'Gmp$'\t'Imk$'\t'Lsk$'\t'Mep$'\t'Mon$'\t'Neu$'\t'Er4$'\t'G1e > rnaHtseqCount_noM_withcoords.mean.txt
tail -n+2 rnaHtseqCount_noM_withcoords.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6, ($7+$8)/2, ($9+$10)/2, ($11+$12)/2, ($13+$14)/2, ($15+$16)/2, ($17+$18)/2, ($19+$20)/2, ($21+$22)/2, $23, ($24+$25)/2, ($26+$27)/2, ($28+$29)/2}' >> rnaHtseqCount_noM_withcoords.mean.txt
python count_to_tpm.py rnaHtseqCount_noM_withcoords.mean.txt rnaTPM_noM_withcoords.mean.txt

#split TPMfile into train, test, and ref inputs saving as npz's
python split_input_rna.py ../../usevision/train_genes_uniq.bed ../../usevision/test_genes_uniq.bed ../../usevision/ref_genes_uniq.bed ../data/RNAseq/rnaTPM_noM_withcoords.mean.txt

#awk'd out chrY

python prepare_expression.py rnaTPM_train_noY.txt 36225 trainTPM.npz
python prepare_expression.py rnaTPM_test_noY.txt 3975 testTPM.npz
python prepare_expression.py rnaTPM_ref_noY.txt 1541 refTPM.npz

#prepare proportion of states in 75kbp window around train, test, and ref TSSs
python window_around_TSS.py rnaTPM_train_noY.txt train_TSS_window.txt 75000
python window_around_TSS.py rnaTPM_test_noY.txt test_TSS_window.txt 75000
python window_around_TSS.py rnaTPM_ref_noY.txt ref_TSS_window.txt 75000

date; time bedtools intersect -a train_TSS_window.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > trainTSS_windows_int_all.bed
date; time bedtools intersect -a train_TSS_window.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > trainTSS_windows_v_all.bed

date; time bedtools intersect -a test_TSS_window.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > testTSS_windows_int_all.bed
date; time bedtools intersect -a test_TSS_window.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > testTSS_windows_v_all.bed

date; time bedtools intersect -a ref_TSS_window.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > refTSS_windows_int_all.bed
date; time bedtools intersect -a ref_TSS_window.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > refTSS_windows_v_all.bed

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

date; time python ccre_prop_in_state.py trainTSS_windows_collapse_{}.bed 37550 trainTSS_window_state_prop.npz
date; time python ccre_prop_in_state.py testTSS_windows_collapse_{}.bed 4218 testTSS_window_state_prop.npz
date; time python ccre_prop_in_state.py refTSS_windows_collapse_{}.bed 1541 refTSS_window_state_prop.npz

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
