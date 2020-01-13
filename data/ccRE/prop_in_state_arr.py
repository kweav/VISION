#!/usr/bin/env python3

import numpy as np

cellTypesOI = ["Lsk", "Hpc7", "Cmp", "Mep", "G1e", "Er4", "Cfue", "Eryad", "Eryfl", "Cfum", "Imk", "Mk", "Gmp", "Mon", "Neu", "Clp", "Nk", "B", "Tcd4", "Tcd8"]
overall_file = 'VISIONmusHem_ccREs_wct_filterkw.bed'
intersect_file = 'VISIONmusHem_ccREs_wct_filterkw_int_{}.bed'

ccRE_index = {}
index_ccRE = {}
for i, line in enumerate(open(overall_file)):
    fields = line.strip('\r\n').split('\t')
    ccRE_index[(fields[0], fields[1], fields[2], fields[3])] = i
    index_ccRE[i] = (fields[0], fields[1], fields[2], fields[3])


num_ccRE = len(ccRE_index.keys())
num_cellsOI = len(cellTypesOI)
cellType_index = dict(zip(cellTypesOI, range(num_cellsOI)))
states = list(range(27))
num_states = len(states)
cellTypeSpecifics = np.zeros((num_ccRE,num_cellsOI,num_states))

def get_overlap(range_1, range_2): #range_1 is the one you want to divide by its length
    range_1_len = range_1[1] - range_1[0]
    overlap_region = [max(range_1[0], range_2[0]), min(range_1[1], range_2[1])]
    len_overlap_region = overlap_region[1] - overlap_region[0]
    overlap = round(len_overlap_region/range_1_len, 5)
    return(overlap)

for cellTypeOI in cellTypesOI:
    ctIndex = cellType_index[cellTypeOI]
    for line in open(intersect_file.format(cellTypeOI)):
        fields=line.strip('\r\n').split('\t')
        this_ccRE = (fields[0], fields[1], fields[2], fields[3])
        this_ccRE_index = ccRE_index[this_ccRE]
        state_range = [int(fields[6]), int(fields[7])]
        ccRE_range = [int(fields[1]), int(fields[2])]
        overlap = get_overlap(ccRE_range, state_range)
        state = int(fields[8])
        cellTypeSpecifics[this_ccRE_index, ctIndex, state] += overlap

fileToWriteTo = open('VISIONmusHem_ccREs_wct_filterkw_prop.bed', 'w+')
for k in range(num_ccRE):
    chr, start, stop, access = index_ccRE[k]
    newline = chr + '\t' + start + '\t' + stop + '\t' + access + '\t'
    subfield = ''
    for cellTypeOI in cellTypesOI:
        subfield += cellTypeOI + '='
        ctIndex = cellType_index[cellTypeOI]
        states_act_props = np.nonzero(cellTypeSpecifics[k,ctIndex])[0]
        act_props = cellTypeSpecifics[k,ctIndex,states_act_props]
        for state, prop in zip(states_act_props, act_props):
            subfield += '(' + str(state) + '_' + str(prop) + ')'
        subfield += ';'
    newline += subfield + '\n'
    fileToWriteTo.write(newline)

# actualProps = np.nonzero(cellTypeSpecifics)
# print(cellTypeSpecifics.shape)
# print(actualProps[0].shape)
# print(actualProps[1].shape)
# print(actualProps[2].shape)
#
# print(len(list(set(actualProps[0]))))
#
# #print(np.nonzero(np.diff(actualProps[0]))[0])
