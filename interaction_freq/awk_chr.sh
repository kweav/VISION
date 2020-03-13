#!/bin/bash

source activate basic

head -n 1 $1 > interaction_freq_all_header.txt

for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX
do
  awk -v var=$CHR '{if ($1 == var) {print}}' $1 > mouse_${CHR}_allInfo.txt
  python plot_int_freq.py $CHR mouse_${CHR}_allInfo.txt interaction_freq_all_header.txt 2048
done

python plot_full_int_freq.py original_paired_{}.npz
python plot_num_paired.py original_paired_{}.npz

conda deactivate
