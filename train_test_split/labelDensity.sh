#!/bin/bash

#SBATCH --job-name=labelDensity
#SBATCH --time=72:0:0
#SBATCH --ntasks-per-node=48
#SBATCH --mem=1000000M
#SBATCH --partition=lrgmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kweave23@jhu.edu
#SBATCH --output=labelDensity.out

export PATH=$PATH:$HOME/miniconda2/bin
source activate /home-3/kweave23@jhu.edu/miniconda2/envs/basic

function bedtools_it {
  bedtools intersect -wa -wb -a $1 -b $2 > $3
  bedtools intersect -v -a $1 -b $2 > $4
}

WINDOWSBED=/home-3/kweave23@jhu.edu/work/users/kweave23/data/windows.bed
PATH1=/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8Seg
PATH2=/home-3/kweave23@jhu.edu/work/users/kweave23/data/

bedtools_it $WINDOWSBED ${PATH1}B.bed ${PATH1}B.window.bed ${PATH1}B.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Cd4.bed ${PATH1}Cd4.window.bed ${PATH1}Cd4.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Cd8.bed ${PATH1}Cd8.window.bed ${PATH1}Cd8.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Cfue.bed ${PATH1}Cfue.window.bed ${PATH1}Cfue.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Cfum.bed ${PATH1}Cfum.window.bed ${PATH1}Cfum.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Clp.bed ${PATH1}Clp.window.bed ${PATH1}Clp.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Cmp.bed ${PATH1}Cmp.window.bed ${PATH1}Cmp.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Er4.bed ${PATH1}Er4.window.bed ${PATH1}Er4.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Eryad.bed ${PATH1}Eryad.window.bed ${PATH1}Eryad.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Eryfl.bed ${PATH1}Eryfl.window.bed ${PATH1}Eryfl.complement.bed
bedtools_it $WINDOWSBED ${PATH1}G1e.bed ${PATH1}G1e.window.bed ${PATH1}G1e.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Gmp.bed ${PATH1}Gmp.window.bed ${PATH1}Gmp.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Hpc7.bed ${PATH1}Hpc7.window.bed ${PATH1}Hpc7.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Imk.bed ${PATH1}Imk.window.bed ${PATH1}Imk.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Lsk.bed ${PATH1}Lsk.window.bed ${PATH1}Lsk.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Mep.bed ${PATH1}Mep.window.bed ${PATH1}Mep.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Mk.bed ${PATH1}Mk.window.bed  ${PATH1}Mk.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Mon.bed ${PATH1}Mon.window.bed ${PATH1}Mon.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Neu.bed ${PATH1}Neu.window.bed ${PATH1}Neu.complement.bed
bedtools_it $WINDOWSBED ${PATH1}Nk.bed ${PATH1}Nk.window.bed ${PATH1}Nk.complement.bed
bedtools_it $WINDOWSBED ${PATH2}HPC7_Rep1_10K_cscore.bg ${PATH2}HPC7_Rep1.window.bed ${PATH2}HPC7_Rep1.complement.bed
bedtools_it $WINDOWSBED ${PATH2}HPC7_Rep2_10K_cscore.bg ${PATH2}HPC7_Rep2.window.bed ${PATH2}HPC7_Rep2.complement.bed
bedtools_it $WINDOWSBED ${PATH2}GSM2514768_10K_cscore.bg ${PATH2}GSM2514768.window.bed ${PATH2}GSM2514768.complement.bed

# python /home-3/kweave23@jhu.edu/work/users/kweave23/train_test_split/labelDensity.py --IDEAS_files \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegB.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegCd4.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegCd8.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegCfue.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegCfum.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegClp.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegCmp.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegEr4.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegEryad.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegEryfl.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegG1e.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegGmp.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegHpc7.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegImk.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegLsk.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegMep.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegMk.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegMon.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegNeu.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegNk.window.bed" --complement_files "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegB.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegCd4.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegCd8.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegCfue.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegCfum.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegClp.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegCmp.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegEr4.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegEryad.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegEryfl.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegG1e.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegGmp.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegHpc7.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegImk.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegLsk.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegMep.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegMk.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegMon.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegNeu.complement.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8SegNk.complement.bed" --LAD_files "/home-3/kweave23@jhu.edu/work/users/kweave23/data/HPC7_Rep1.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/HPC7_Rep2.window.bed" \
# "/home-3/kweave23@jhu.edu/work/users/kweave23/data/GSM2514768.window.bed" > printStatements.out

conda deactivate
