#!/bin/bash

source activate bwa_env

FILE_START=SRR52958
FILE_END_1F=_1.fastq
FILE_END_2F=_2.fastq
FILE_END_1B=_1_bwa.bam
FILE_END_2B=_2_bwa.bam
F=first
S=sec

for NUM in $(seq $1 $2)
do
  date; time bwa mem -t 20 -T 30 -P bwa_index_ /project/vision/HiC_Modeling/Fastq/$FILE_START$NUM$FILE_END_1F | samtools view -b -o $FILE_START$NUM$FILE_END_1B
  date; time bwa mem -t 20 -T 30 -P bwa_index_ /project/vision/HiC_Modeling/Fastq/$FILE_START$NUM$FILE_END_2F | samtools view -b -o $FILE_START$NUM$FILE_END_2B
  date; time samtools view -f 0x40 -b -o $FILE_START$NUM$F$FILE_END_1B $FILE_START$NUM$FILE_END_1B
  date; time samtools view -f 0x40 -b -o $FILE_START$NUM$F$FILE_END_2B $FILE_START$NUM$FILE_END_2B
  date; time samtools view -f 0x80 -b -o $FILE_START$NUM$S$FILE_END_1B $FILE_START$NUM$FILE_END_1B
  date; time samtools view -f 0x80 -b -o $FILE_START$NUM$S$FILE_END_2B $FILE_START$NUM$FILE_END_2B
  echo $NUM complete
done

conda deactivate
