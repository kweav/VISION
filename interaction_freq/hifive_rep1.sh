#!/bin/bash

FILE_START=SRR52958
FILE_END_B=_bwa.bam
F=_1
S=_2

source activate hifive_env
date; time hifive hic-complete express -B mm10_DpnII_frag.bed --re DpnII --genome mm10 -i 650 -f 3 -j 500 -n 75 \
-S ${FILE_START}58$F$FILE_END_B ${FILE_START}58$S$FILE_END_B \
-S ${FILE_START}59$F$FILE_END_B ${FILE_START}59$S$FILE_END_B \
-S ${FILE_START}60$F$FILE_END_B ${FILE_START}60$S$FILE_END_B \
-S ${FILE_START}61$F$FILE_END_B ${FILE_START}61$S$FILE_END_B \
-S ${FILE_START}62$F$FILE_END_B ${FILE_START}62$S$FILE_END_B \
-P G1e_rep1_
conda deactivate
