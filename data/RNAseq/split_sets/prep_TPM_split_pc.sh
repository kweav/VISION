#!/bin/bash

#prepare TPM files from count data
awk '{if ($1 != "chrM") {print}}' rnaHtseqCountsall_withcoordinates.0.txt > rnaHtseqCount_noM_withcoords.txt

echo chr$'\t'start$'\t'end$'\t'gene$'\t'gene_type$'\t'strand$'\t'Cfue$'\t'Cfum$'\t'Cmp$'\t'Eryfl$'\t'Gmp$'\t'Imk$'\t'Lsk$'\t'Mep$'\t'Mon$'\t'Neu$'\t'Er4$'\t'G1e > rnaHtseqCount_noM_withcoords.mean.txt
tail -n+2 rnaHtseqCount_noM_withcoords.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6, ($7+$8)/2, ($9+$10)/2, ($11+$12)/2, ($13+$14)/2, ($15+$16)/2, ($17+$18)/2, ($19+$20)/2, ($21+$22)/2, $23, ($24+$25)/2, ($26+$27)/2, ($28+$29)/2}' >> rnaHtseqCount_noM_withcoords.mean.txt
python count_to_tpm.py rnaHtseqCount_noM_withcoords.mean.txt rnaTPM_noM_withcoords.mean.txt

#split TPMfile into train, test, and ref inputs saving as npz's
python split_input_rna.py ~/taylorLab/usevision/train_genes_uniq.bed ~/taylorLab/usevision/test_genes_uniq.bed ~/taylorLab/usevision/ref_genes_uniq.bed rnaTPM_noM_withcoords.mean.txt

#awk out chrY
awk '{if ($1 != "chrY") {print}}' rnaTPM_train.txt > rnaTPM_train_noY.txt
awk '{if ($1 != "chrY") {print}}' rnaTPM_test.txt > rnaTPM_test_noY.txt
awk '{if ($1 != "chrY") {print}}' rnaTPM_ref.txt > rnaTPM_ref_noY.txt

#awk out just protein coding genes
head -n 1 rnaTPM_train_noY.txt > rnaTPM_train_pc.txt
awk '{if ($5 == "protein_coding") {print}}' rnaTPM_train_noY.txt >> rnaTPM_train_pc.txt
head -n 1 rnaTPM_test_noY.txt > rnaTPM_test_pc.txt
awk '{if ($5 == "protein_coding") {print}}' rnaTPM_test_noY.txt >> rnaTPM_test_pc.txt
head -n 1 rnaTPM_ref_noY.txt > rnaTPM_ref_pc.txt
awk '{if ($5 == "protein_coding") {print}}' rnaTPM_ref_noY.txt >> rnaTPM_ref_pc.txt

#make npz files
./prepare_expression.py rnaTPM_train_pc.txt 19012 trainTPM_pc.npz
./prepare_expression.py rnaTPM_test_pc.txt 1967 testTPM_pc.npz
./prepare_expression.py rnaTPM_ref_pc.txt 857 refTPM_pc.npz
