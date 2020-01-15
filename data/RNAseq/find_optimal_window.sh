#!/bin/bash

#SBATCH --job-name=TSS_WINDOW5
#SBATCH --time=1:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kweave23@jhu.edu
#SBATCH --output=outfile20200115_TSS_window.out

source activate /work-zfs/jtayl139/kweave23/my_envs/basic

IDEAS=/home-3/kweave23@jhu.edu/data/kweave23/IDEAS/ideasVisionV20p8Seg
BED=.bed

#tail -n +2 $1 | awk '{print $1, $2, $2}' OFS='\t' > TSSs.bed
for WIN in 57500 60000 62500 65000 70000 75000
do
  for CT in B Cd4 Cd8 Cfue Cfum Clp Cmp Er4 Eryad Eryfl G1e Gmp Hpc7 Imk Mk Lsk Mep Mon Neu Nk
  do
    bedtools window -a TSSs.bed -b $IDEAS$CT$BED -w $WIN -v > TSS_windows_v_${WIN}_$CT$BED
  done
done

wc -l TSS_windows_v_*.bed

conda deactivate
