#!/bin/bash

#PATH1=/Users/kateweaver/taylorLab/VISION/data/IDEAS/onlyTranscribed/ideasVisionV20p8Seg*.transcribed.bed
#WINDOWS=/Users/kateweaver/taylorLab/VISION/data/IDEAS/windows.bed
#PATH2=/Users/kateweaver/taylorLab/VISION/data/IDEAS/only2/ideasVisionV20p8Seg*.bed.o2.bed

#PATH1P1=/project/vision/Data/IDEAS/ideasVisionV20p8Seg
#PATH1P2=bed
PATH2P1=/Users/kateweaver/taylorLab/VISION/data/IDEAS/only2/ideasVisionV20p8Seg
PATH2P2=bed.o2.bed
#WINDOWS1=/project/vision/Data/IDEAS/windows_100K_10K.bed
#WINDOWS2=/project/vision/Data/IDEAS/windows_75000_5000.bed
#WINDOWS3=/project/vision/Data/IDEAS/windows_55000_5000.bed
WINDOWS3=/Users/kateweaver/taylorLab/VISION/data/IDEAS/windows_bed_files/windows_55000_5000.bed
#WINDOWTYPE1=100K_10K
#WINDOWTYPE2=75K_5K
WINDOWTYPE3=55K_5K

function bedtools_it {
  bedtools intersect -wa -wb -a $1 -b $2 > $3.$5.$4.bed
  bedtools intersect -v -a $1 -b $2 > $3.$6.$4.bed
}

for TYPE in B Cd4 Cd8 Cfue Cfum Clp Cmp Er4 Eryad Eryfl G1e Gmp Hpc7 Imk Lsk Mep Mk Mon Neu Nk
do
  #FILEOI=$PATH1P1$TYPE.$PATH1P2
  #APPEND1_1=window
  #APPEND1_2=complement
  FILEOI2=$PATH2P1$TYPE.$PATH2P2
  APPEND2_1=o2.window
  APPEND2_2=o2.complement
  echo $FILEOI2
  #bedtools_it $WINDOWS1 $FILEOI $PATH1P1$TYPE $WINDOWTYPE1 $APPEND1_1 $APPEND1_2
  #bedtools_it $WINDOWS2 $FILEOI $PATH1P1$TYPE $WINDOWTYPE2 $APPEND1_1 $APPEND1_2
  #bedtools_it $WINDOWS3 $FILEOI $PATH1P1$TYPE $WINDOWTYPE3 $APPEND1_1 $APPEND1_2
  bedtools_it $WINDOWS3 $FILEOI2 $PATH2P1$TYPE $WINDOWTYPE3 $APPEND2_1 $APPEND2_2
done

#for FILE1 in $PATH1
#do
#  bedtools_it $WINDOWS $FILE1
#done

#for FILE2 in $PATH2
#do
#  bedtools_it $WINDOWS $FILE2
#done
