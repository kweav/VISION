#!/usr/bin/env python3

import numpy as np
import argparse as ap

parser = ap.ArgumentParser(description='Looking at the Label Representation in windows, specifically amount heterochromatin and related amount Transcribed')
parser.add_argument('--window', action='store', nargs=1, type=int, required=False, default=[75000], help='size of sliding window')
parser.add_argument('--slide', action='store', nargs=1, type=int, required=False, default=[5000], help='size of slide')
args=parser.parse_args()
window = args.window[0]
slide = args.slide[0]

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                'chr18', 'chr19', 'chrX', 'chrY']

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

windowsFile=open('windows_{}_{}.bed'.format(window, slide), 'w+')
for chromosome in chromosomes:
    chrEnd = chrSizes[chromosome]
    finalWindowStart = chrEnd - window
    windowArray = np.arange(0,finalWindowStart, slide)
    for value in windowArray:
        windowsFile.write(chromosome + '\t' + str(value) + '\t' + str(value+window) +'\n')
    windowsFile.write(chromosome + '\t' + str(finalWindowStart) + '\t' + str(chrEnd) + '\n')
windowsFile.close()
