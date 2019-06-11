#!/bin/bash

#SBATCH --job-name=cscore
#SBATCH --time=72:0:0
#SBATCH --ntasks-per-node=24
#SBATCH --partition=shared
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kweave23@jhu.edu
#SBATCH --output=cscore.out

export PATH=$PATH:$HOME/miniconda2/bin
source activate /home-3/kweave23@jhu.edu/miniconda2/envs/basic

python /home-3/kweave23@jhu.edu/work/users/kweave23/train_test_split/cscoreDict.py --LAD_files "/home-3/kweave23@jhu.edu/work/users/kweave23/data/HPC7_Rep1.window.bed" \
"/home-3/kweave23@jhu.edu/work/users/kweave23/data/HPC7_Rep2.window.bed" \
"/home-3/kweave23@jhu.edu/work/users/kweave23/data/GSM2514768.window.bed" > printStatements_cscore.out

conda deactivate
