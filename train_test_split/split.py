#!/usr/bin/env python3

import argparse as ap
import logging
import datetime
import numpy as np
import sklearn.model_selection as sms
import matplotlib.pyplot as plt
import collections
from statistics import mean

logging.basicConfig(filename='split.log',level=logging.DEBUG)
logging.info(datetime.datetime.now())

parser = ap.ArgumentParser(description='split the pruned windows into a training and test set stratifying on GC content')
parser.add_argument('--prunedWindows', action='store', nargs=1, type=str, required=True, help='file in bed format with gc content in final field')
parser.add_argument('--decimals', action='store', nargs=1, type=int, required=False, default=[2], help='number of decimal places to round GC content to')
parser.add_argument('--test_size', action='store', nargs=1, type=float, required=False, default=[0.1], help='proportion for test size')
args=parser.parse_args()
prunedWindows = args.prunedWindows[0]
decimals = args.decimals[0]
test_size = args.test_size[0]

'''make list to be split and make list to stratify split with'''
windows = []
gc = []
for line in open(prunedWindows):
    fields=line.strip('\r\n').split('\t')
    windows.append((fields[0], fields[1], fields[2], fields[-1]))
    gc.append(float(fields[-1]))
logging.info("file values added to list " + str(datetime.datetime.now()))

gc=np.array(gc)
gc=np.round(gc, decimals=decimals)
logging.info("rounded gc values to {} decimals ".format(str(decimals)) + str(datetime.datetime.now()))

'''split'''
windowsTrain, windowsTest = sms.train_test_split(windows, test_size = test_size, random_state = 47, stratify = gc)
logging.info("split windows stratifying on gc with a test_size of {} ".format(test_size) + str(datetime.datetime.now()))

'''record split and assess split for....'''
trainLens = []
trainGCs = []
trainChrs = []
f=open("train.bed", "w+")
for i, value in enumerate(windowsTrain, 1):
    chr, start, stop, gcV = value
    f.write(chr + '\t' + start + '\t' + stop + '\n')
    trainLens.append(int(stop)-int(start))
    trainGCs.append(float("{0:.2f}".format(float(gcV))))
    trainChrs.append(chr)
numTrainTotal = i
f.close()
logging.info("wrote values to train.bed")

testLens = []
testGCs = []
testChrs = []
f = open("test.bed", 'w+')
for i, value in enumerate(windowsTest, 1):
    chr, start, stop, gcV = value
    f.write(chr + '\t' + start + '\t' + stop + '\n')
    testLens.append(int(stop)-int(start))
    testGCs.append(float("{0:.2f}".format(float(gcV))))
    testChrs.append(chr)
numTestTotal = i
f.close()
logging.info("wrote values to test.bed")

''' ....length/amount of genome'''
fig, (ax1, ax2) = plt.subplots(1,2)
ax1.hist(trainLens, bins=25)
ax1.set_xlabel("Length of Region")
ax1.set_title("Train Set")
ax1.xaxis.set_tick_params(rotation=45)
ax2.hist(testLens, bins=25)
ax2.set_xlabel("Length of Region")
ax2.set_title("Test Set")
ax2.xaxis.set_tick_params(rotation=45)
plt.tight_layout()
fig.savefig("RegionLens.png")
plt.close(fig)
logging.info("RegionLens.png plotted " + str(datetime.datetime.now()))

meanTrainLen = mean(trainLens)
meanTestLen = mean(testLens)

print("Mean Training Region Length: ", meanTrainLen)
logging.info("mean training region length " + str(meanTrainLen))
print("Mean Test Region Length: ", meanTestLen)
logging.info("mean test region length " + str(meanTestLen))

numUniqTrain = len(set(trainGCs))
TrainChr = collections.Counter(trainChrs)
sumTrainLens = sum(trainLens)
print("sum of lens in train set: ", sumTrainLens)
logging.info("sum of lens in train set " + str(sumTrainLens))

numUniqTest = len(set(testGCs))
TestChr = collections.Counter(testChrs)
sumTestLens = sum(testLens)
print("sum of lens in test set: ", sumTestLens)
logging.info("sum of lens in test set " + str(sumTestLens))

print("true length proportion: ", sumTestLens/sumTrainLens)
logging.info("true length proportion " + str(sumTestLens/sumTrainLens))

'''....GC content represented '''
fig, (ax1, ax2) = plt.subplots(1,2)
ax1.hist(trainGCs, bins=numUniqTrain)
ax1.set_xlabel("GC content")
ax1.set_title("Train set")
ax2.hist(testGCs, bins=numUniqTest)
ax2.set_xlabel("GC content")
ax2.set_title("Test set")
fig.savefig("GCsplit.png")
plt.close(fig)
logging.info("GCsplit.png plotted " + str(datetime.datetime.now()))

'''... chromosomes represented '''
chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                'chr18', 'chr19', 'chrX', 'chrY']

trainHeight = []
testHeight = []
for chromosome in chromosomes:
    trainHeight.append(TrainChr[chromosome]/numTrainTotal)
    testHeight.append(TestChr[chromosome]/numTestTotal)
fig, ax = plt.subplots()
ax.bar(np.arange(len(chromosomes))-0.2, trainHeight, width=0.4)
ax.bar(np.arange(len(chromosomes))+0.2, testHeight, width=0.4)
ax.set_xticks(np.arange(len(chromosomes)))
ax.set_xticklabels(chromosomes, rotation = 90)
plt.legend(labels=["Train", "Test"])
fig.savefig("Chrsplit.png")
plt.close(fig)
logging.info("Chrsplit.png plotted " + str(datetime.datetime.now()))
