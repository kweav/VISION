# Script order & Brief description of believed purpose

## 1) Training 

**Output**: `Rdat` files for each chromosome {1,2,...,19, X}, each gene group {1,2,3,4}, and each lessone cell type {1,2,...,12} 

             vision_rna_tss2k_ccreunit.{chr}.{gene_group}.{lessone}.Rdat
A) `rungenecre.R`
  - Driver for `genecre.R`
  - loops through chromosomes, gene groups, and lessones
  - saves Rdat files
  
B) `genecre.R`
  - subselection procedure
  - returns dictionary with pairings, selections, correlations
  
## 2) Extract Pairing information and Reproducibility of Pairing among lessones
## 3) Calculate adjR2 values
## 4) Calculate final Betas
## 5) Get eRP tracks
