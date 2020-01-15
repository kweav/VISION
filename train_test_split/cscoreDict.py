#!/usr/bin/env python3

import argparse as ap
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import logging
import datetime
import pickle

logging.basicConfig(filename='cScoreDict.log',level=logging.DEBUG)

parser = ap.ArgumentParser(description='Consider IDEAS label density and correlation with LAD score')
parser.add_argument('--LAD_files', action='store', nargs='+', type=str, required=True, help='cscore/LAD_score files')
parser.add_argument('--window', action='store', nargs=1, type=int, required=False, default=[100000], help='size of sliding window')
parser.add_argument('--slide', action='store', nargs=1, type=int, required=False, default=[10000], help='size of slide')
args=parser.parse_args()
LAD_files = args.LAD_files
window = args.window[0]
slide = args.slide[0]

windowsLocLAD ={}
for file in LAD_files:
    logging.info('beginning {}'.format(file))
    logging.info(datetime.datetime.now())
    for line in open(file):
        fields=line.strip('\r\n').split('\t')
        chr=fields[0]
        if chr not in windowsLocLAD:
            windowsLocLAD[chr]={}
        window_start = int(fields[1])
        window_end = int(fields[2])
        if (window_start, window_end) not in windowsLocLAD[chr]:
            windowsLocLAD[chr][(window_start, window_end)] = []
        windowsLocLAD[chr][(window_start, window_end)].append(float(fields[6]))
    logging.info('finished {}'.format(file))
    logging.info(datetime.datetime.now())

pickle_out = open("windowsLocLAD_Dict.pickle", "wb")
pickle.dump(windowsLocLAD, pickle_out)
pickle_out.close()
