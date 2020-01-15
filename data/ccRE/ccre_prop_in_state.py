#!/usr/bin/env python3

import numpy as np
import sys

'''Usage: ./ccre_prop_in_state.py input_file_base num_locs outfile'''
'''Usage: ./ccre_prop_in_state.py VISIONmusHem_ccRE_collapse_{}.bed 205019 ccre_state_prop.npz'''
'''Usage: ./ccre_prep_in_state.py TSS_windows_collapse_{}.bed 41814 TSS_window_state_prop.npz'''

'''
Tab-delimited 4 column input file structure; 20 different files hense an input_file_base:
chr locOI_start locOI_stop  comma_separated_list_of_overlapping_segments_from_bedtools_groupby_collapse
This script aims to find the proportion of each IDEAS state for each given genomic location of interest
Actions taken:
1) If only one segment overlaps a locOI:
    1.1) The proportion overlap of the state is added for that locOI and cellType
2) If more than one segment overlaps a locOI, states and proportion of overlap are stored in ordered numpy arrays for each overlapping segment unless that segment is 10kbp+ and state 0:
    2.1) If the rounded to 2 decimal places sum of the overlap proportions is less than or equal to 1.0, the proportion overlap of the each state is added for that locOI and cellType
    2.2) If the rounded to 2 deimcal places sum sum of the overlap proprotions is greater than 1.0, the unique states which are represented is found:
        2.2.1) If a single unique state is represented, the actual covered proportion is set for that single state
        2.2.2) If more than one unique state is represented, an ad hoc standardization technique is used such that the final sum of proportions is 1
Finally save 3 arrays:
1 array has all the proportions
1 array has the cell type indexing
1 array has the ccRE indexing
'''



def get_overlap(range_1, range_2, ccRE_len):
    range_1_len = range_1[1] - range_1[0]
    range_2_len = range_2[1] - range_2[0]
    overlap_region = [max(range_1[0], range_2[0]), min(range_1[1], range_2[1])]
    len_overlap_region = overlap_region[1] - overlap_region[0]
    overlap = len_overlap_region/ccRE_len
    return (overlap)


cellTypes = ['Lsk', 'Hpc7', 'Cmp', 'Mep', 'G1e', 'Er4', 'Cfue', 'Eryad', 'Eryfl',
            'Cfum', 'Imk', 'Mk', 'Gmp', 'Mon', 'Neu', 'Clp', 'Nk', 'B', 'Cd4', 'Cd8']
cellN = len(cellTypes)
cellType_index = dict(zip(cellTypes, range(cellN)))
states = np.arange(27)
stateN = states.shape[0]
ccREn = int(sys.argv[2]) #205019

file_base = sys.argv[1] #'VISIONmusHem_ccRE_collapse_{}.bed'

props = np.zeros((ccREn, cellN, stateN))

cre_index = {}
for k, ct in enumerate(cellTypes):
    for i, line in enumerate(open(file_base.format(ct))):
        fields=line.strip('\r\n').split('\t')
        chr, start, stop = fields[0:3]
        len_ccRE = int(stop) - int(start) + 1
        overlapping = fields[3].split(',')
        if len(overlapping) == 1: #just a single overlapping IDEAS state. Add it.
            istart, istop, state = overlapping[0].split('_')
            overlap = get_overlap([int(start), int(stop)], [int(istart), int(istop)], len_ccRE)
            if k == 0: #add based on i since ccRE index dictionary based off of first cell type
                props[i, cellType_index[ct], int(state)] += overlap
            else: #add based on definite ccRE index from index dictionary
                props[cre_index[(chr, start, stop)], cellType_index[ct], int(state)] += overlap
        else:
            overlaps = np.zeros(len(overlapping), dtype=np.float32)
            states_represented = np.full(len(overlapping), -1, dtype=np.int32)
            cre_locs = np.zeros((len_ccRE, len(overlapping)))
            for j, instance in enumerate(overlapping):
                istart, istop, state = instance.split('_')
                ilen = int(istop) - int(istart) + 1
                if ilen >= 10000 and int(state) == 0: #this happens 196458 times for 181258 lines, and is reached in all 20 cell types
                    continue
                else:
                    states_represented[j] = int(state)
                    overlaps[j] = get_overlap([int(start), int(stop)], [int(istart), int(istop)], len_ccRE)
                    cre_locs[max(0, int(istart)-int(start)): min(len_ccRE-1, int(istop)-int(start)+1), j] = 1
            rounded_sum = np.around(np.sum(overlaps),2) #can round before looking at valid_locs because the overlap will be zero
            valid_locs = states_represented != -1
            if rounded_sum <= 1.0: #in the case that the total overlaps are less than 1, add it. If total overlaps are one, add it. If above one, don't add it yet
                for value_o, value_s in zip(overlaps[valid_locs], states_represented[valid_locs]):
                    if k == 0: #add based on i since ccRE index dictionary based off of first cell type
                        props[i, cellType_index[ct], value_s] += value_o
                    else: #add based on definite ccRE index from index dictionary
                        props[cre_index[(chr, start, stop)], cellType_index[ct], value_s] += value_o
            else: #want to check what in overlapping is overlapping each other (a single state or a mixture of states?)
                unique_states = np.unique(states_represented[valid_locs])
                if unique_states.shape[0] == 1:
                    actually_covered_prop = (1 - np.sum(np.sum(cre_locs, axis=1)==0)/cre_locs.shape[0]) #finding locations that don't have any coverage and subtracting that proporition from 1
                    if k == 0: #add based on i since ccRE index dictionary based off of first cell type
                        props[i, cellType_index[ct], unique_states[0]] += actually_covered_prop
                    else: #add based on definite ccRE index from index dictionary
                        props[cre_index[(chr, start, stop)], cellType_index[ct], unique_states[0]] += actually_covered_prop

                else:
                    overlaps[valid_locs] /= np.sum(overlaps[valid_locs]) #ad hoc standardization so that the sum of proportions is 1; this will penalize overlapping.
                    for value_o, value_s in zip(overlaps[valid_locs], states_represented[valid_locs]):
                        if k == 0: #add based on i since ccRE index dictionary based off of first cell type
                            props[i, cellType_index[ct], value_s] += value_o
                        else: #add based on definite ccRE index from index dictionary
                            props[cre_index[(chr, start, stop)], cellType_index[ct], value_s] += value_o

        if k == 0: #build index dictionary for ccREs
            cre_index[(chr, start, stop)] = i

print(np.sum(props > 1.01))
print(props[props > 1.01])

cellIndex = np.empty(cellN, dtype=np.object)
for ct in cellTypes:
    cellIndex[cellType_index[ct]] = ct

ccREIndex = np.empty((ccREn, 3), dtype=np.object)
for key in cre_index:
    ccREIndex[cre_index[key], :] = key

outfile = sys.argv[3]
f = open(outfile, 'wb')
np.savez(f, props = props, cellIndex = cellIndex, ccREIndex = ccREIndex)
f.close()
