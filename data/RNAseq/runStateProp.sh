#!/bin/bash

#SBATCH --job-name=TSS_initial
#SBATCH --time=12:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kweave23@jhu.edu
#SBATCH --output=outfile20200115_TSS_initial.out

'''Usage: sbatch runStateProp.sh ~/data/kweave23/RNAseq/rnaTPM.txt 75000 41814'''

source activate /work-zfs/jtayl139/kweave23/my_envs/basic

tail -n +2 $1 | awk '{print $1, $2}' OFS='\t' > TSSs.txt
date; time python window_around_TSS.py TSSs.txt $2

date; time bedtools intersect -a TSSs_withWindow.bed -b /home-3/kweave23@jhu.edu/data/kweave23/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > TSS_windows_int_all.bed
date; time bedtools intersect -a TSSs_withWindow.bed -b /home-3/kweave23@jhu.edu/data/kweave23/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > TSS_windows_v_all.bed

IDEAS=/home-3/kweave23@jhu.edu/data/kweave23/IDEAS/ideasVisionV20p8Seg
BED=.bed
for CT in B Cd4 Cd8 Cfue Cfum Clp Cmp Er4 Eryad Eryfl G1e Gmp Hpc7 Imk Mk Lsk Mep Mon Neu Nk
do
  date; time awk -v var=$IDEAS$CT$BED '{if ($4 == var) {print}}' TSS_windows_int_all.bed > TSS_windows_int_$CT$BED
  date; time awk '{print $1, $2, $3, $6"_"$7"_"$8}' OFS='\t' TSS_windows_int_$CT$BED > TSS_windows_istate_$CT$BED
  date; time bedtools groupby -i TSS_windows_istate_$CT$BED -g 1,2,3 -c 4 -o collapse > TSS_windows_collapse_$CT$BED
done

date; time python ccre_prop_in_state.py TSS_windows_collapse_{}.bed $3 TSS_window_state_prop.npz

conda deactivate
