#!/bin/bash

awk '{if ($1 != "chrM") {print}}' rnaHtseqCountsall_withcoordinates.0.txt > rnaHtseqCount_noM_withcoords.txt

echo chr$'\t'start$'\t'end$'\t'gene$'\t'gene_type$'\t'strand$'\t'Cfue$'\t'Cfum$'\t'Cmp$'\t'Eryfl$'\t'Gmp$'\t'Imk$'\t'Lsk$'\t'Mep$'\t'Mon$'\t'Neu$'\t'Er4$'\t'G1e > rnaHtseqCount_noM_withcoords.mean.txt
tail -n+2 rnaHtseqCount_noM_withcoords.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6, ($7+$8)/2, ($9+$10)/2, ($11+$12)/2, ($13+$14)/2, ($15+$16)/2, ($17+$18)/2, ($19+$20)/2, ($21+$22)/2, $23, ($24+$25)/2, ($26+$27)/2, ($28+$29)/2}' >> rnaHtseqCount_noM_withcoords.mean.txt
python count_to_tpm.py rnaHtseqCount_noM_withcoords.mean.txt rnaTPM_noM_withcoords.mean.txt

echo chr$'\t'start$'\t'end$'\t'gene$'\t'gene_type$'\t'strand$'\t'Cfue$'\t'Cfum$'\t'Cmp$'\t'Eryfl$'\t'Gmp$'\t'Imk$'\t'Lsk$'\t'Mep$'\t'Mon$'\t'Neu$'\t'Er4$'\t'G1e > rnaHtseqCount_noM_withcoords.sum.txt
tail -n+2 rnaHtseqCount_noM_withcoords.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6, ($7+$8), ($9+$10), ($11+$12), ($13+$14), ($15+$16), ($17+$18), ($19+$20), ($21+$22), $23, ($24+$25), ($26+$27), ($28+$29)}' >> rnaHtseqCount_noM_withcoords.sum.txt
python count_to_tpm.py rnaHtseqCount_noM_withcoords.sum.txt rnaTPM_noM_withcoords.sum.txt

python split_input_rna.py ../../usevision/train_genes_uniq.bed ../../usevision/test_genes_uniq.bed ../../usevision/ref_genes_uniq.bed ../data/RNAseq/rnaTPM_noM_withcoords.mean.txt
