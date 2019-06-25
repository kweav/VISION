#!/bin/bash
SEP=_
END=.txt
COC=coverage
TAPREF=T
TOP=1000
WINDOW=55000
SLIDE=5000
for AFFIRM in no yes
do
  for NUM in 15 18
  do
    for RATIO in 1.35 1.5 2.0
    do
      #NEW=/home/kweave23/VISION_train_test_split/labelRepresentation/75K_5K_{$COC}_{$AFFIRM}_{$TAPREF}_{$NUM}_{$RATIO}_{$TOP}
      NEW=/home/kweave23/VISION_train_test_split/labelRepresentation/55K_5K_$COC$SEP$AFFIRM$SEP$TAPREF$SEP$NUM$SEP$RATIO$SEP$TOP
      mkdir $NEW
      cd $NEW
      python ./../labelRepresentation.py --window $WINDOW --slide $SLIDE --COC $COC --countQ $AFFIRM --TApref $TAPREF --numCellTypes $NUM --thresholdRatio $RATIO --thresholdTop $TOP
      awk '{print $1}' potentialWindows_params_55K_5K_$COC$SEP$AFFIRM$SEP$TAPREF$SEP$NUM$SEP$RATIO$SEP$TOP$END | sort | uniq -c > chrom.txt
      echo finished with {$COC}_{$AFFIRM}_{$TAPREF}_{$NUM}_{$RATIO}
    done
  done
done
