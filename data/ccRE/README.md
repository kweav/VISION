# ccRE and accessibility data

## VISIONmusHem_ccREs_filterkw.bed
5 field tab delimited file with chr, start, stop,then accessibility by celltype and epigentic/IDEAS state by celltype
Formatted as such:
chr1    4769915 4770171 0_0_1_0_0_0_0_1_1_1_0_0_0_1_1_0_0_0     25;13;13;5;13;13;13;13;13;13;25;13;7;3;13;13;7;25;25;24;

>./format_ccRE.py VISIONmusHem_ccREs.bed > VISIONmusHem_ccREs_filterkw.bed

## VISIONmusHem_ccREs_wct_filterkw.bed
5 field tab delimited file with chr, start, stop,then accessibility by celltype and epigentic/IDEAS state by celltype. Differs from `VISIONmusHem_ccREs_filterkw.bed` by explicitly specifying cellType. Formatted as such:
chr1    3071299 3071451 Lsk=0;Hpc7=1;Cmp=0;Mep=0;G1e=0;Er4=0;Cfue=0;Eryad=0;Eryfl=0;Cfum=0;Imk=0;Gmp=0;Mon=0;Neu=0;Nk=0;B=0;Tcd4=0;Tcd8=0;      Lsk=0;Hpc7=9;Cmp=0;Mep=0;G1e=0;Er4=0;Cfue=0;Eryad=0;Eryfl=0;Cfum=0;Imk=0;Mk=0;Gmp=0;Mon=0;Neu=0;Clp=0;Nk=0;B=0;Tcd4=0;Tcd8=0;

>./format_ccRE_wct.py VISIONmusHem_ccREs.bed > VISIONmusHem_ccREs_wct_filterkw.bed

## VISIONmusHem_ccREs_filterkw2.bed
3 field tab delimited file with the first field being a unique identifier for genomic location, 2nd field being accessibility by celltype and 3rd field being epigenetic/IDEAS state by celltype.
Formatted as such:
Ah..48c87bh..48c97b     0_0_1_0_0_0_0_1_1_1_0_0_0_1_1_0_0_0     25;13;13;5;13;13;13;13;13;13;25;13;7;3;13;13;7;25;25;24;

> ./genomeLoc_to_ID.py --input_file VISIONmusHem_ccREs_filterkw.bed --ccRE_chr_field 0 --ccRE_start_field 1 --ccRE_end_field 2 --output_file VISIONmusHem_ccREs_filterkw2_1.bed
> cat VISIONmusHem_ccREs_filterkw2_1.bed | cut -f6,7,8 > VISIONmusHem_ccREs_filterkw2.bed

## VISIONmusHem_ccREs_wct_filterkw2.bed
3 field tab delimited file with the first field being a unique identifier for genomic location, 2nd field being accessibility by celltype and 3rd field being epigenetic/IDEAS state by celltype. Differs from `VISIONmusHem_ccREs_filterkw2.bed` by explicitly specifying cellType. Formatted as such:
Ah..2edd43h..2edddb     Lsk=0;Hpc7=1;Cmp=0;Mep=0;G1e=0;Er4=0;Cfue=0;Eryad=0;Eryfl=0;Cfum=0;Imk=0;Gmp=0;Mon=0;Neu=0;Nk=0;B=0;Tcd4=0;Tcd8=0;      Lsk=0;Hpc7=9;Cmp=0;Mep=0;G1e=0;Er4=0;Cfue=0;Eryad=0;Eryfl=0;Cfum=0;Imk=0;Mk=0;Gmp=0;Mon=0;Neu=0;Clp=0;Nk=0;B=0;Tcd4=0;Tcd8=0;

> ./genomeLoc_to_ID.py --input_file VISIONmusHem_ccREs_wct_filterkw.bed --ccRE_chr_field 0 --ccRE_start_field 1 --ccRE_end_field 2 --output_file VISIONmusHem_ccREs_wct_filterkw2_1.bed
> cat VISIONmusHem_ccREs_wct_filterkw2_1.bed | cut -f6,7,8 > VISIONmusHem_ccREs_wct_filterkw2.bed

### more on unique identifier for genomic locations
The script `genomeLoc_to_ID.py` produces 19 character unique IDs where
    -   the first character points to chr (dictionary `chrEncode`)
    -   the rest of the ID can be split on 'h.' to get the start and ending genome locations in a pseudo hexadecimal form; used 'h.' to separate start and end since h alone could be true hexadecimal
        -   First, replace all filler '.' with ''; these fillers were used so that all IDs were the same length
        -   Second, add '0x' to the front to produce the true hexadecimal
        -   Third, convert hexadecimal back to int
        -   Decoding script is `ID_to_genomeLoc.py`

>chrEncode = {'chr1':'A',
            'chr2':'B',
            'chr3':'C',
            'chr4':'D',
            'chr5':'E',
            'chr6':'F',
            'chr7':'G',
            'chr8':'H',
            'chr9':'I',
            'chr10':'J',
            'chr11':'K',
            'chr12':'L',
            'chr13':'M',
            'chr14':'N',
            'chr15':'O',
            'chr16':'P',
            'chr17':'Q',
            'chr18':'R',
            'chr19':'S',
            'chrX':'T',
            'chrY':'U'}

### ccREs within windows of genes (use `bedtools window -w`)
**tl;dr:**
  - `genes_ccRE_{}_window.bed`:
  **tab delimited file with gene location, gene ID, TPM, ccRE location, ccRE accessibility by cellType, ccRE epigenetic IDEAS state by cell type where the ccRE is within {} bp of the gene**
  chr3    108107280       108146146       ENSMUSG00000000001.4    Lsk=81.94;Cmp=65.34;Gmp=96.05;Mep=72.91;Cfue=8.61;Eryad=16.54;Cfum=103.60;Imk=94.43;Mon=75.31;Neu=59.91;G1e=240.12;Er4=130.20;  chr3    108112723       108112988
         Lsk=0;Hpc7=0;Cmp=0;Mep=1;G1e=0;Er4=0;Cfue=0;Eryad=0;Eryfl=0;Cfum=1;Imk=0;Gmp=0;Mon=0;Neu=0;Nk=0;B=0;Tcd4=0;Tcd8=0;      Lsk=1;Hpc7=1;Cmp=1;Mep=9;G1e=1;Er4=8;Cfue=1;Eryad=1;Eryfl=1;Cfum=9;Imk=1;Mk=1;Gmp=1;Mon=6;Neu=1;Clp=1;Nk=1;B=1;Tcd4=1;Tcd8=1;
  - `genes_ccRE_ID_{}_window.bed`:
  **tab delimited file with gene location, gene ID, TPM, ccRE unique identifier, ccRE accessibility by cellType, ccRE epigenetic IDEAS state by cell type where the ccRE is within {} bp of the gene**
    chr3    108107280       108146146       ENSMUSG00000000001.4    Lsk=81.94;Cmp=65.34;Gmp=96.05;Mep=72.91;Cfue=8.61;Eryad=16.54;Cfum=103.60;Imk=94.43;Mon=75.31;Neu=59.91;G1e=240.12;Er4=130.20;  Ch.671aef2h.671afdb     Lsk=0;Hpc7=0;Cmp=0;Mep=0;G1e=0;Er4=0;Cfue=0;Eryad=0;Eryfl=0;Cfum=0;Imk=0;Gmp=0;Mon=1;Neu=0;Nk=0;B=0;Tcd4=0;Tcd8=0;      Lsk=1;Hpc7=1;Cmp=1;Mep=1;G1e=1;Er4=8;Cfue=1;Eryad=1;Eryfl=1;Cfum=1;Imk=1;Mk=1;Gmp=6;Mon=9;Neu=4;Clp=1;Nk=1;B=1;Tcd4=1;Tcd8=1;
  - 3 windows used:
    - within 200bp
    - within 1kb 
    - within 1Mb (but excluding the within 1kb)



- Used bedtools window to find ccREs within 200bp of each gene, 1kb, and 1Mb
  - `bedtools window -w 200 -a scriptseq3.v3.kw2.IDlocexp.bed -b ../ccRE/VISIONmusHem_ccREs_filterkw.bed > genes_ccRE_200bp_window.bed`
  - `bedtools window -w 1000 -a scriptseq3.v3.kw2.IDlocexp.bed -b ../ccRE/VISIONmusHem_ccREs_filterkw.bed > genes_ccRE_1kb_window.bed`
  - `bedtools window -w 1000000 -a scriptseq3.v3.kw2.IDlocexp.bed -b ../ccRE/VISIONmusHem_ccREs_filterkw.bed > genes_ccRE_1Mb_window.bed`  
- But really I only wanted the 200bp and 1kb, and then the 1Mb difference 1kb.
  - `sort -k1,1 -k2,2 -k3,3 -k6,6 -k7,7 -k8,8 genes_ccRE_1kb_window.bed > genes_ccRE_1kb_window_sorted.bed`
  - `sort -k1,1 -k2,2 -k3,3 -k6,6 -k7,7 -k8,8 genes_ccRE_1Mb_window.bed > genes_ccRE_1Mb_window_sorted.bed`
  - `cat genes_ccRE_1Mb_window_sorted.bed >> genes_ccRE_1Mb_1kb.bed`
  - `cat genes_ccRE_1kb_window_sorted.bed >> genes_ccRE_1Mb_1kb.bed` 
  - `sort -k1,1 -k2,2 -k3,3 -k6,6 -k7,7 -k8,8 -k4,4 -k9,9 -k10,10 genes_ccRE_1Mb_1kb.bed | uniq -c > genes_ccRE_1Mb_1kb_suc.bed`
  - `awk '{if ($1 == 1) {print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}}' genes_ccRE_1Mb_1kb_suc.bed > genes_ccRE_1Mb_diff_1kb.bed`
  - `wc -l genes_ccRE_1*.bed`
         >152976 genes_ccRE_1kb_window.bed
         >152976 genes_ccRE_1kb_window_sorted.bed
        >8998627 genes_ccRE_1Mb_1kb.bed
        >8845651 genes_ccRE_1Mb_1kb_suc.bed
        >8692675 genes_ccRE_1Mb_diff_1kb.bed
        >8845651 genes_ccRE_1Mb_window.bed
        >8845651 genes_ccRE_1Mb_window_sorted.bed

    - We can see here that 1Mb "by itself" is 8845651 lines (really this by itself is just with a single count of 1kb included)
    - The concatenated file is 1Mb ("by itself") plus 1kb by itself (8845651 + 152976 = 8998627 ) or in other words 2 copies of the 1kb
    - The suc file which has been sorted and piped to uniq -c is back to 8845651 or what you would suspect as it has collapsed the 2 1kb entries into single 1kb entries again
    - Finally after 'awk'ing out anything with a count of 2, the file is 8692675 lines or 8845651-152976 suggesting that all 1kb entries have been succesfully removed and the 'diff' file has >1kb, <1Mb
  - `./genomeLoc_to_ID.py --input_file genes_ccRE_1Mb_diff_1kb.bed --output_file genes_ccRE_ID_1Mb_diff_1kb.bed`
  - `./genomeLoc_to_ID.py --input_file genes_ccRE_1kb_window.bed    --output_file genes_ccRE_ID_1kb_window.bed`                                                  
  - `./genomeLoc_to_ID.py --input_file genes_ccRE_200bp_window.bed --output_file genes_ccRE_ID_200bp_window.bed`
