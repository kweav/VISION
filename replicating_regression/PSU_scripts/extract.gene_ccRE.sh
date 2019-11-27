module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools

working_dir='/storage/home/gzx103/group/projects/vision/rna'
filename='vision_rna_tss2k_ccreunit'

cd $working_dir

time Rscript get_gene_groups.R $filename

### extract information
for i in {1..12}
do
	echo $i
	cat $filename'.chr'*'.'*'.'$i'.gene_ccRE.txt' > 'vision_rna.gene_ccRE.c'$i'.selected.all.txt'
done

### get reproducibility in leave-one-out
time python check_reproducibility.py

###
time cat gencode.vM4.annotation.gtf | awk -F '"' -v OFS='\t' '{print $1,$2,$10}' | awk -F '\t' -v OFS='\t' '{if ($7=="+") print $1"_"$4, $11":"$10; else if ($7=="-") print $1"_"$5, $11":"$10}' > gencode.vM4.tss.tss2gene0.txt
###

### get gene name
time python get_gene_nameid.py


sort -k3,3 vision_rna.gene_ccRE.selected.all.reprod_count.WithName1209.txt | cut -f1,2,3,4,5,6,8,10 > vision_rna.gene_ccRE.selected.all.reprod_count.WithNameSorted.txt



