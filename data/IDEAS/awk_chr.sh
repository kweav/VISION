#!/bin/bash

PATH1P1=/project/vision/Data/IDEAS/intersect_window_100K_10K/ideasVisionV20p8Seg
PATH1P2=.window.100K_10K.bed

PATH2P1=/project/vision/Data/IDEAS/complement_window_100K_10K/ideasVisionV20p8Seg
PATH2P2=.complement.100K_10K.bed

PATH3P1=/project/vision/Data/IDEAS/intersect_window_75000_5000/ideasVisionV20p8Seg
PATH3P2=.window.75K_5K.bed

PATH4P1=/project/vision/Data/IDEAS/complement_window_75000_5000/ideasVisionV20p8Seg
PATH4P2=.complement.75K_5K.bed

function awk_it {
  for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY
  do
    awk '{if ($1 == "'$CHR'") {print}}' $1 > $2.$CHR.$3
  done
}

function call_it {
  for TYPE in B Cd4 Cd8 Cfue Cfum Clp Cmp Er4 Eryad Eryfl G1e Gmp Hpc7 Imk Lsk Mep Mk Mon Neu Nk
  do
    FILEOI=$1$TYPE$2
    awk_it $FILEOI $1$TYPE $2
  done
}

call_it $PATH1P1 $PATH1P2
call_it $PATH2P1 $PATH2P2
call_it $PATH3P1 $PATH3P2
call_it $PATH4P1 $PATH4P2
