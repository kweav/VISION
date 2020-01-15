#!/usr/bin/env python3

'''Usage ./window_around_TSS.py TSSs.txt window_size'''

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
            'chrY':91744698,
            'chr18':90702639,
            'chr19':61431566}

import sys
window = int(sys.argv[2])
fileToWriteTo = open('TSSs_withWindow.bed', 'w+')


for line in open(sys.argv[1]):
    chr, TSS = line.strip('\r\n').split('\t')
    lineToWrite = '{}\t{}\t{}\n'.format(chr, max(1, round(float(TSS))-window), min(chrSizes[chr], round(float(TSS))+window))
    fileToWriteTo.write(lineToWrite)

fileToWriteTo.close()
