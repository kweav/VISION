#!/bin/bash

for window in 70000 65000 60000 55000 50000 45000 40000 35000 30000 25000 20000 15000 10000 9500 9000 8500 8000 7500 7000 6500 6000 5500 5000
do
  python window_around_TSS.py rnaTPM_train_pc.txt train_TSS_${window}_pc.txt $window
  python window_around_TSS.py rnaTPM_test_pc.txt test_TSS_${window}_pc.txt $window
  python window_around_TSS.py rnaTPM_ref_pc.txt ref_TSS_${window}_pc.txt $window
  echo $window
  echo train
  date; time bedtools intersect -a train_TSS_${window}_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames | wc -l
  echo test
  date; time bedtools intersect -a test_TSS_${window}_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames | wc -l
  echo ref
  date; time bedtools intersect -a ref_TSS_${window}_pc.txt -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames | wc -l
done
