#!/bin/bash

FILE_START=SRR52958
FILE_END_B=_bwa.bam
F=first
S=sec

conda activate bwa_env
for NUM in 73 74 75 76 77
do
  date; time samtools view -f 0x40 -b -o $FILE_START$NUM$F$FILE_END_B $FILE_START$NUM$FILE_END_B
  date; time samtools view -f 0x80 -b -o $FILE_START$NUM$S$FILE_END_B $FILE_START$NUM$FILE_END_B
  echo $NUM complete
done
conda deactivate

conda deactivate hifive_env
date; time hifive hic-complete express -B mm10_DpnII_frag.bed --re DpnII --genome mm10 -d -i 0 \
-S $FILE_START{73}$F$FILE_END_B $FILE_START{73}$S$FILE_END_B \
-S $FILE_START{74}$F$FILE_END_B $FILE_START{74}$S$FILE_END_B \
-S $FILE_START{75}$F$FILE_END_B $FILE_START{75}$S$FILE_END_B \
-S $FILE_START{76}$F$FILE_END_B $FILE_START{76}$S$FILE_END_B \
-S $FILE_START{77}$F$FILE_END_B $FILE_START{77}$S$FILE_END_B \
-P G1e_rep2_
conda deactivate
