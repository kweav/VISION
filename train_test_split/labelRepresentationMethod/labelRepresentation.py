#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import logging
import datetime
import argparse as ap

logging.basicConfig(filename='labelRepresentation.log',level=logging.DEBUG)
logging.info('began job' + str(datetime.datetime.now()))

parser = ap.ArgumentParser(description='Looking at the Label Representation in windows, specifically amount heterochromatin and related amount Transcribed')
parser.add_argument('--window', action='store', nargs=1, type=int, required=False, default=[75000], help='size of sliding window')
parser.add_argument('--slide', action='store', nargs=1, type=int, required=False, default=[5000], help='size of slide')
parser.add_argument('--COC', action='store', nargs=1, type=str, required=False, default=['coverage'], help='coverage or count; if count, will look at number of labels; if coverage, will look at amount of window covered by label' )
parser.add_argument('--countQ', action='store', nargs=1, type=str, required=False, default=['yes'], help='yes or no; if yes, will include quiescent in consideration; if no, will not include quiescent in considerations')
parser.add_argument('--TApref', action='store', nargs=1, type=str, required=False, default=['TA'], help='T or TA for limiting just transcribed or limiting transcribed & active')
parser.add_argument('--numCellTypes', action='store', nargs=1, type=int, required=False, default=[18], help='minimum number of cell types to have threshold values to make it a candidate region')
parser.add_argument('--thresholdRatio', action='store', nargs=1, type=float, required=False, default=[3.0], help='minimum threshold ratio of difference like such: amount_covered_Heterochromatin/total_covered_window - amount_covered_Transcribed/total_covered_window = diff_H-T/total_covered_window ~ H-T=D. If D>0, thresholdRatio*T < H')
parser.add_argument('--thresholdTop', action='store', nargs=1, type=int, required=False, default=[7500], help='after argsort of maximum dead windows and then finding the highest ratios, I take a slice of atleast this many to compare among cell types for that chromosome')
args=parser.parse_args()
windowG = args.window[0]
windowAb = str(windowG).replace("0", "") + "K"
if windowAb == '1K':
    windowAb = '100K'
slideG = args.slide[0]
slideAb = str(slideG).replace("0", "") + "K"
if slideAb == '1K':
    slideAb = '10K'
coc = args.COC[0] #count or coverage
countQ = args.countQ[0] #no or yes for including quiescent labels in overall count/coverage
TApref = args.TApref[0] #T or TA
numCellTypes_threshold = args.numCellTypes[0]
thresholdRatio = args.thresholdRatio[0]
thresholdTop = args.thresholdTop[0]

logging.info('window: ' + str(windowG) + '\nslide: ' + str(slideG) + '\nCOC: ' + coc + '\ncountQ: ' + countQ + '\nTApref: ' + TApref + '\nnumCellTypes: ' + str(numCellTypes_threshold) + '\n_thresholdRatio_*T < H: ' + str(thresholdRatio) + '\nthresholdTop: ' + str(thresholdTop))

#ideas_window = '/home-3/kweave23@jhu.edu/work/users/kweave23/data/ideasVisionV20p8Seg{}.window.bed.{}.bed' #cell type, chromosome
ideas_window = '/project/vision/Data/IDEAS/intersect_window_{}_{}/chrAwk/ideasVisionV20p8Seg{}.{}..window.{}_{}.bed' #str(windowG), str(slideG), cellType, chromosome, windowAb, slideAb
#ideas_window = '/Users/kateweaver/taylorLab/VISION/data/IDEAS/chr_awk/ideasVisionV20p8Seg{}.{}.75K_5K.window.head25.bed' #cellTypes, chromosome

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                 'chr18', 'chr19', 'chrX', 'chrY']

cellTypes = ['B', 'Cd4', 'Cd8', 'Cfue', 'Cfum', 'Clp', 'Cmp', 'Er4', 'Eryad', 'Eryfl', 'G1e', 'Gmp', 'Hpc7', 'Imk', 'Lsk', 'Mep', 'Mk', 'Mon', 'Neu', 'Nk']

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

#send labels to 1 of 8 specific overall colors/categories
labelDict = { 2 : "Heterochromatin",
              1 : "Transcribed",
              8 : "Transcribed",
              14 : "Transcribed",
              17 : "Transcribed",
              25 : "Transcribed",
              4 : "Enhancer",
              5 : "Enhancer",
              3 : "Polycomb/Bivalent",
              16 : "Polycomb/Bivalent",
              20 : "Polycomb/Bivalent",
              22 : "Polycomb/Bivalent",
              10 : "CRM",
              11 : "CRM",
              15 : "CRM",
              18 : "CRM",
              21 : "CRM",
              24 : "CRM" ,
              6 : "Active",
              12 : "Active",
              19 : "Active",
              23 : "Active",
              7 : "CTCF/NucleaseAccessible",
              9 : "CTCF/NucleaseAccessible",
              13 : "CTCF/NucleaseAccessible",
              26 : "CTCF/NucleaseAccessible",
              0 : "Quiescent"}

#translate label type to depth index
depthIndex = {'Heterochromatin':0,
              'Transcribed':1,
              'Enhancer':2,
              'CRM': 3,
              'Active':4,
              'Polycomb/Bivalent':5,
              'CTCF/NucleaseAccessible':6,
              'Quiescent':7}

labelTypes = list(set(labelDict.values()))
if countQ == 'no':
    labelTypes.remove('Quiescent')

'''make 3 index dictionaries
1) num_windows provides the number of windows (value) per chromosome (key)
2) window_index provides the index (value) of a chromosome and (window_start, window_end) tuple (key and subkey). This index correponds to 2nd dimension location in later array
3) index_window provides the (window_start, window_end) tuple (value) for a given chromosome and and index (key and subkey)'''
num_windows = {}
window_index = {}
index_window = {}
for chromosome in chromosomes:
    num_windows[chromosome] = 0
    chrEnd = chrSizes[chromosome]
    finalWindowStart = chrEnd - windowG
    windowArray = np.arange(0,finalWindowStart, slideG)
    window_index[chromosome] = {}
    index_window[chromosome] = {}
    for i, value in enumerate(windowArray):
        window_index[chromosome][(value, value+windowG)] = i
        index_window[chromosome][i] = (value, value+windowG)
        num_windows[chromosome] += 1
    window_index[chromosome][(finalWindowStart, chrEnd)] = i+1
    index_window[chromosome][i+1] = (finalWindowStart, chrEnd)
    num_windows[chromosome] += 1
max_num_windows = max(num_windows.values())

logging.info('index dictionaries initialized ' + str(datetime.datetime.now()))

'''make a storage dictionary:
1) windowOI - store how many times a window for a given chromosome (subkey and key) meets the criteria for a cell type (value is an int ranging from 1-20)'''
windowsOI = {}
for chromosome in chromosomes:
    windowsOI[chromosome] = {}

logging.info('storage dictionaries intialized ' + str(datetime.datetime.now()))

'''fill an array (CHR_array) with count/coverage of label types
fill an array (window_coverage) with actual bps covered by labels in the window'''
potentialWindows = open('potentialWindows_params_{}_{}_{}_{}_{}_{}_{}_{}.txt'.format(str(windowG), str(slideG), coc, countQ, TApref, str(numCellTypes_threshold), str(thresholdRatio), str(thresholdTop)), 'w+')
potentialWindows.write('Chromosome\tPotential_CHB_Start\tPotential_CHB_End\tAverage(Heterochromatin_slice_allCellTypes_thatWindow)\tStdev(H_slice)\tAverage(Transcribed_slice_allCellTypes_thatWindow)\tStdev(T_slice)\n') #remove this line later using tail -n +2
for j, chromosome in enumerate(chromosomes):
    CHR_array = np.full((len(cellTypes), num_windows[chromosome],len(labelTypes)), 0)
    window_coverage = np.full((len(cellTypes), num_windows[chromosome],1), 0)
    for k, cellType in enumerate(cellTypes):
        ideas_window_file = open(ideas_window.format(str(windowG), str(slideG), cellType, chromosome, windowAb, slideAb))
        for line in ideas_window_file:
            fields= line.strip('\r\n').split('\t')
            window_start = int(fields[1])
            window_end = int(fields[2])
            label_start = int(fields[4])
            label_end = int(fields[5])
            amount_covered = label_end - label_start
            column_index = window_index[chromosome][(window_start, window_end)]
            label = int(fields[6])
            if countQ == 'no' and label == 0:
                pass
            else:
                if coc == 'count':
                    CHR_array[k,column_index, depthIndex[labelDict[label]]] += 1
                    window_coverage[k,column_index,0] += 1
                elif coc == 'coverage':
                    CHR_array[k,column_index, depthIndex[labelDict[label]]] += amount_covered
                    window_coverage[k, column_index,0] += amount_covered
    logging.info('finished parsing {} for {} cellType into arrays'.format(chromosome, cellType) + str(datetime.datetime.now()))

    '''store information in storage dictionary'''
    fraction_covered = np.where(window_coverage > 0, np.divide(CHR_array, window_coverage), 0)
    logging.info('normalized counts/coverage by total count or coverage in window ' + str(datetime.datetime.now()))
    for k, cellType in enumerate(cellTypes):
        maximum_dead_indices = np.argsort(fraction_covered[k, :, depthIndex['Heterochromatin']])[::-1] #sort high to low
        H = fraction_covered[k,:,depthIndex['Heterochromatin']][maximum_dead_indices]
        if TApref == 'T':
            T = corresponding_T = fraction_covered[k,:,depthIndex['Transcribed']][maximum_dead_indices]
        elif TApref == 'TA':
            T = fraction_covered[k,:,depthIndex['Transcribed']][maximum_dead_indices] + fraction_covered[k,:,depthIndex['Active']][maximum_dead_indices]
        D = H - T
        indices_diff = D > 0
        H_imp = H[indices_diff]
        T_imp = T[indices_diff]
        mdi_imp = maximum_dead_indices[indices_diff]
        ratioT_imp = T_imp * thresholdRatio
        indices_ratio = H_imp > ratioT_imp
        true_max_dead_indices = mdi_imp[indices_ratio]
        logging.info(str(true_max_dead_indices.shape[0]) + " for " + cellType)
        if true_max_dead_indices.shape[0] > thresholdTop: #because I sorted high to low earlier, grabbing regions that meet D and thresholdRatio criteria but also have higher levels of H to begin with
            cut_true_max_dead_indices=true_max_dead_indices[:thresholdTop]
        else:
            cut_true_max_dead_indices = true_max_dead_indices
        for value in cut_true_max_dead_indices:
            window = index_window[chromosome][value]
            if window not in windowsOI[chromosome]:
                windowsOI[chromosome][window] = 0
            windowsOI[chromosome][window] += 1
    logging.info('finished storing how many cell types a window on {} meets my criteria '.format(chromosome) + str(datetime.datetime.now()))
    '''now go through and find potential areas according to the numCellTypes_threshold and print these potential CHB areas to a file and plot where on chromosome the areas are''' #will merge overlapping regions later
    potWindowsOI=[]
    fig, ax = plt.subplots()
    plt.suptitle("{} Heterochromatic Regions OI".format(chromosome))
    ax.axvline(0, linestyle='--')
    ax.axvline(1, linestyle='--')
    ax.set_xlabel('Relative chromosome location')
    ax.set_ylabel('Region Number')
    for i in windowsOI[chromosome]:
        actual_numCellTypes = windowsOI[chromosome][i]
        if actual_numCellTypes >= numCellTypes_threshold:
            potWindowsOI.append((i, actual_numCellTypes))
    for i, value in enumerate(potWindowsOI):
        potWindowOI, actual_numCellTypes = value
        potWindow_start, potWindow_end = potWindowOI
        ax.scatter(int(potWindow_start)/chrSizes[chromosome], i, marker='>')
        ax.scatter(int(potWindow_end)/chrSizes[chromosome], i, marker='<')
        window_index_val = window_index[chromosome][potWindowOI]
        H_slice = fraction_covered[:,window_index_val,0]
        T_slice = fraction_covered[:,window_index_val,1]
        potentialWindows.write(chromosome + '\t' + str(potWindow_start) + '\t' + str(potWindow_end) + '\t' + str(np.average(H_slice)) + '\t' + str(np.std(H_slice)) + '\t' + str(np.average(T_slice)) + '\t' + str(np.std(T_slice)) +'\n')
    plotName = 'where_on_{}_with_params_{}_{}_{}_{}_{}_{}_{}_{}.png'.format(chromosome, str(windowG), str(slideG), coc, countQ, TApref, str(numCellTypes_threshold), str(thresholdRatio), str(thresholdTop))
    fig.savefig(plotName)
    plt.close(fig)
    logging.info('plotted location of regions of interest on {}. The plots is saved as {}'.format(chromosome, plotName) + str(datetime.datetime.now()))
    logging.info('finished assessment of regions of interest that met my criteria on {} '.format(chromosome) + str(datetime.datetime.now()))
