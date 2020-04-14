#!/bin/bash

for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX
do
  python use_clustered_regress.py --where comp --chroms $CHR --n_thresh 30 --dist_thresh 1000 &
done
