#!/usr/bin/env python3

import argparse as ap
import numpy as np

parser = ap.ArgumentParser(description='order the mean GC content')
parser.add_argument('--tenMB_GC', action='store', nargs=1, type=str, required = True, help='file from computeMeanGC.py with avg 10MB GC content')
args=parser.parse_args()
tenMB_GC = args.tenMB_GC[0]

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

locations = []
gcContents = []

previousEnd = 0
for i, line in enumerate(open(tenMB_GC)):
    fields=line.strip('\r\n').split(' ')
    if fields[-1] != 'final':
        stop = int(fields[2])+1000000
        previousEnd = stop
        start = stop - 10000000
        location = (fields[0], start, stop)
    else:
        stop = chrSizes[fields[0]]
        location = (fields[0], previousEnd, stop)
    locations.append(location)
    gcMean = float(fields[1])
    gcContents.append(gcMean)

locations = np.array(locations)
gcContents = np.array(gcContents)
wouldSort = np.argsort(gcContents)
sorted_locations = locations[wouldSort]
sorted_gcContents = gcContents[wouldSort]

fileToWriteTo = open('sorted_10MB_GC.txt', 'w+')
for loc, gc in zip(sorted_locations, sorted_gcContents):
    chr, start, stop = loc
    fileToWriteTo.write(chr + '\t' + str(start) + '\t' + str(stop) + '\t' + str(gc) + '\n')
