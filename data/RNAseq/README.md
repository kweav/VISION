## RNAseq data files
`filter_TPM_files.py` was used to filter `amitData/{cellType}.amit.genes.tsv` to produce `amitData/{cellType}.amit.geneID_TPM.tab` and `amitData/{cellType}.amit.geneID.txt`
a similar/but simplified version of the same script (`filter_TPM_files_gen.py`) was used to filter `scriptseq3.v3.tab` to `scriptseq3.v3.kw1.tab`


### scriptseq3.v3.kw2.tab
tab delimited file with 4 fields (specifically) chr, start, stop, and TPM of 12 Hematopoietic cell types with the final field being formatted as such: Lsk=81.94;Cmp=65.34;Gmp=96.05;Mep=72.91;Cfue=8.61;Eryad=16.54;Cfum=103.60;Imk=94.43;Mon=75.31;Neu=59.91;G1e=240.12;Er4=130.20;

>./geneID_to_loc.py --gencode ~/mm10_genome/gencode.vM4.annotation.gtf --RNAseq scriptseq3.v3.kw1.tab


### scriptseq3.v3.kw2.IDlocexp.bed
tab delimited file with 5 fields (specifically) chr, start, stop, geneID, and TPM of 12 Hematopoietic cell types with the final field being formatted as such:
Lsk=81.94;Cmp=65.34;Gmp=96.05;Mep=72.91;Cfue=8.61;Eryad=16.54;Cfum=103.60;Imk=94.43;Mon=75.31;Neu=59.91;G1e=240.12;Er4=130.20;

>./geneID_loc_exp.py --gencode ~/mm10_genome/gencode.vM4.annotation.gtf --RNAseq scriptseq3.v3.kw1.tab

### amit.cellTypes_withloc.tab
tab delimited file with 4 fields (specifically) chr, start, stop, and TPM of 4 adaptive immunity cell types with the final field being formatted as such:
B=8.91;Nk=17.10;Tcd4=11.89;Tcd8=8.23;

>./geneID_to_loc.py --gencode ~/mm10_genome/gencode.vM4.annotation.gtf --RNAseq amitData/\*.amit.geneID_TPM.tab

### amit.cellTypes_IDlocexp.bed
tab delimited file with 5 fields (specifically) chr, start, stop, geneID, and TPM of 4 adaptive immunity cell types with the final field being formatted as such:
B=8.91;Nk=17.10;Tcd4=11.89;Tcd8=8.23;

>./geneID_loc_exp.py --gencode ~/mm10_genome/gencode.vM4.annotation.gtf --RNAseq amitData/\*.amit.geneID_TPM.tab

### both_RNAseq.tab
tab delimited file with 4 fields (specifically) chr, start, stop, and TPM of 16 Hematopoietic cell types (including adaptive immunity) with the final field being formatted such:
Lsk=140.06;Cmp=53.73;Gmp=37.27;Mep=54.14;Cfue=6.31;Eryad=10.93;Cfum=27.87;Imk=112.75;Mon=218.05;Neu=231.16;G1e=49.25;Er4=84.56;B=13.16;Nk=47.41;Tcd4=18.21;Tcd8=16.65

>./geneID_to_loc2.py --gencode ~/mm10_genome/gencode.vM4.annotation.gtf --RNAseq scriptseq3.v3.kw1.tab amitData/*.amit.geneID_TPM.tab
