#!/usr/bin/env python

import argparse as ap
import numpy as np
import subprocess
import logging
import datetime

logging.basicConfig(filename='removeReference.log',level=logging.DEBUG)

logging.info(datetime.datetime.now())

def findLeaps(arrayOI, window_start, window_end, interval_start, interval_end):
    intervalLeaps = []
    if interval_start != window_start:
        LeapStart = interval_start
    else:
        LeapStart = window_start
    for i, value in enumerate(arrayOI):
        if i !=0:
            previousValue = arrayOI[i-1]
            if value != previousValue+1:
                LeapEnd = previousValue
                intervalLeaps.append((LeapStart, LeapEnd))
                LeapStart= value
            if value == interval_end:
                intervalLeaps.append((LeapStart, interval_end))

    return(intervalLeaps)

parser = ap.ArgumentParser(description='remove blacklist regions from sorted 10MB GC')
parser.add_argument('--region_GC', action='store', nargs=1, type=str, required = True, help='file with sorted avg GC content for regions') #gcMean_proposed_splits_sorted.bed
parser.add_argument('--blacklist', action='store', nargs=1, type=str, required = True, help='bed file with blacklist regions or reference gene Loci') #referenceLoci_0.5Mbp.txt
args=parser.parse_args()
region_GC = args.region_GC[0]
blacklist = args.blacklist[0]

#use bedtools intersect
bedtools_intersect_out = subprocess.check_output("bedtools intersect -a {} -b {} -wa -wb".format(region_GC, blacklist), shell=True).decode("utf-8").splitlines()
logging.info('intersection complete')
logging.info(datetime.datetime.now())
bedtools_dont_intersect_out = subprocess.check_output("bedtools intersect -a {} -b {} -v".format(region_GC, blacklist), shell=True).decode("utf-8").splitlines()
logging.info('intersection complement complete')
logging.info(datetime.datetime.now())
#check the output from bedtools intersect and then change the start and stop to reflect the removal
#fileToWriteTo2 = open("TEST_bedtools_intersect.txt", 'w+')

fileToWriteTo = open('notSorted_10MB_GC_removedBlacklist_reference.txt', 'w+')
windowsDict = {}
for i, value in enumerate(bedtools_intersect_out): #chr0 window_start1 window_end2 gc3 togetherness4 slides5 seqLen6 chr7 blacklist_start8 blacklist_end9 reference_loci_name10
    #fileToWriteTo2.write(value + '\n')
    fields = value.split("\t")
    chr = fields[0]
    gc = fields[3]
    together = fields[4]

    blacklist_start = int(fields[8])
    blacklist_end = int(fields[9])+1
    window_start = int(fields[1])
    window_end = int(fields[2])+1

    valuesOI = (chr, window_start, window_end, gc, together)
    if valuesOI not in windowsDict:
        windowsDict[valuesOI] = []
        window_array = np.arange(window_start, window_end)
        windowsDict[valuesOI] = window_array.copy()

    blacklist_array = np.arange(blacklist_start, blacklist_end)
    window_array = windowsDict[valuesOI].copy()
    difference_array = np.setdiff1d(window_array, blacklist_array)
    windowsDict[valuesOI] = difference_array.copy()

logging.info('setdiff1d complete')
logging.info(datetime.datetime.now())

for value in bedtools_dont_intersect_out:
    fileToWriteTo.write(value + '\n')

logging.info('recording intersection complement complete')
logging.info(datetime.datetime.now())

for key in windowsDict:
    chr, window_start, window_end, gc, together = key
    window_end += -1
    intervalStart = windowsDict[key][0]
    intervalEnd = windowsDict[key][-1]
    leaps = findLeaps(windowsDict[key], window_start, window_end, intervalStart, intervalEnd)
    for value in leaps:
        start, end = value
        fileToWriteTo.write(chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(gc) + '\t'+together+'\n')

logging.info('recording difference intervals complete')
logging.info(datetime.datetime.now())
