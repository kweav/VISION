#!/bin/bash

FILE_START=SRR52958
FILE_END_S=_bwa.sam
FILE_END_B=_bwa.bam
FOR=for
REV=rev

for NUM in 58 59 60 61 62 73 74 75 76 77
do
  date; time samtools view -f 0x10 -b -o $FILE_START$NUM$REV$FILE_END_B $FILE_START$NUM$FILE_END_B
  date; time samtools view -F 0x10 -b -o $FILE_START$NUM$FOR$FILE_END_B $FILE_START$NUM$FILE_END_B
  date; time samtools view -f 0x10 -h -o $FILE_START$NUM$REV$FILE_END_S $FILE_START$NUM$FILE_END_B
  date; time samtools view -F 0x10 -h -o $FILE_START$NUM$FOR$FILE_END_S $FILE_START$NUM$FILE_END_B
  echo $NUM complete
done
