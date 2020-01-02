#!/bin/bash

FILE_START=SRR52958
FILE_END_B=_bwa.bam
F=_1
S=_2

source activate hifive_env
date; time hifive hic-complete express -B mm10_DpnII_frag.bed --re DpnII --genome mm10 -i 650 -f 3 -j 500 -n 75 \
-S ${FILE_START}73$F$FILE_END_B ${FILE_START}73$S$FILE_END_B \
-S ${FILE_START}74$F$FILE_END_B ${FILE_START}74$S$FILE_END_B \
-S ${FILE_START}75$F$FILE_END_B ${FILE_START}75$S$FILE_END_B \
-S ${FILE_START}76$F$FILE_END_B ${FILE_START}76$S$FILE_END_B \
-S ${FILE_START}77$F$FILE_END_B ${FILE_START}77$S$FILE_END_B \
-P G1e_rep2_
conda deactivate
