#!/bin/bash

FILE_START=SRR52958
FILE_END_B=_bwa.bam
F=first
S=sec

for NUM in 58 59 60 61 62
do
  date; time samtools view -f 0x40 -b -o $FILE_START$NUM$F$FILE_END_B $FILE_START$NUM$FILE_END_B
  date; time samtools view -f 0x80 -b -o $FILE_START$NUM$S$FILE_END_B $FILE_START$NUM$FILE_END_B
  echo $NUM complete
done

date; time hifive hic-complete express -B mm10_DpnII_frag.bed --re DpnII --genome mm10 -d -i 0 \
-S $FILE_START{58}$F$FILE_END_B $FILE_START{58}$S$FILE_END_B \
-S $FILE_START{59}$F$FILE_END_B $FILE_START{59}$S$FILE_END_B \
-S $FILE_START{60}$F$FILE_END_B $FILE_START{60}$S$FILE_END_B \
-S $FILE_START{61}$F$FILE_END_B $FILE_START{61}$S$FILE_END_B \
-S $FILE_START{62}$F$FILE_END_B $FILE_START{62}$S$FILE_END_B \
-P G1e_rep1_
