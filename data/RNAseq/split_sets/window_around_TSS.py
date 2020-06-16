#!/usr/bin/env python3

'''Usage ./window_around_TSS.py input_file output_file window_size'''

chrSizes = {'chr1':195471971,
            'chr2':182113224,
            'chrX':171031299,
            'chr3':160039680,
            'chr4':156508116,
            'chr5':151834684,
            'chr6':149736546,
            'chr7':145441459,
            'chr10':130694993,
            'chr8':129401213,
            'chr14':124902244,
            'chr9':124595110,
            'chr11':122082543,
            'chr13':120421639,
            'chr12':120129022,
            'chr15':104043685,
            'chr16':98207768,
            'chr17':94987271,
            'chr18':90702639,
            'chr19':61431566}

import sys
window = int(sys.argv[3])
output_file = sys.argv[2]
fileToWriteTo = open(output_file, 'w+')


for i,line in enumerate(open(sys.argv[1])):
    if i == 0:
        continue
    fields = line.strip('\r\n').split('\t')
    chr = fields[0]
    if fields[5] == '+':
        TSS = int(fields[1])
    elif fields[5] == '-':
        TSS = int(fields[2])
    lineToWrite = '{}\t{}\t{}\t{}\n'.format(chr, max(1, TSS-window), min(chrSizes[chr], TSS+window), fields[3]) #including the gene name because several have the same TSS but are still different genes
    fileToWriteTo.write(lineToWrite)

fileToWriteTo.close()
