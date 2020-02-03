#!/usr/bin/env python3

'''Usage: ./split_input_rna.py ../../usevision/train_genes_uniq.bed ../../usevision/test_genes_uniq.bed ../../usevision/ref_genes_uniq.bed ../data/RNAseq/rnaTPM_noM_withcoords.mean.txt'''

import sys
train_file, test_file, ref_file, input_file = sys.argv[1:]

chroms = []
for i in range(1,20):
    chroms.append('chr{}'.format(i))
chroms.append('chrX')
chroms.append('chrY')

chrom_gene = {'train':{},
             'test':{},
             'ref':{}}
for chrom in chroms:
    for key in chrom_gene:
        chrom_gene[key][chrom] = set()

for file, key in zip([train_file, test_file, ref_file], ['train', 'test', 'ref']):
    for line in open(file):
        fields = line.strip('\r\n').split('\t')
        chr, gene_name = fields[0], fields[3]
        chrom_gene[key][chr].add(gene_name)

train_input = open('rnaTPM_train.txt', 'a+')
test_input = open('rnaTPM_test.txt', 'a+')
ref_input = open('rnaTPM_ref.txt', 'a+')
num_unaccounted = 0
for i,line in enumerate(open(input_file)):
    if i == 0:
        train_input.write('{}\n'.format(line.strip('\r\n')))
        test_input.write('{}\n'.format(line.strip('\r\n')))
        ref_input.write('{}\n'.format(line.strip('\r\n')))
        continue
    fields = line.strip('\r\n').split()
    chr, gene_name = fields[0], fields[3]
    if gene_name in chrom_gene['train'][chr]:
        train_input.write('{}\n'.format(line.strip('\r\n')))
    elif gene_name in chrom_gene['test'][chr]:
        test_input.write('{}\n'.format(line.strip('\r\n')))
    elif gene_name in chrom_gene['ref'][chr]:
        ref_input.write('{}\n'.format(line.strip('\r\n')))
    else:
        num_unaccounted += 1

train_input.close()
test_input.close()
ref_input.close()
print(num_unaccounted)
