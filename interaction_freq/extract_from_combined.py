#!/usr/bin/env python

import sys
import hifive
import matplotlib.pyplot as plt
import numpy as np

font = {'size'   : 18}
plt.rc('font', **font)

'''Usage: ./extract_from_combined.py all_paired_file interaction_length
          ./extract_from_combined.py ~/taylorLab/VISION/data/eRP_scores_PSU_regression/mouse_chr19_allInfo.txt 2000'''

hic_1 = hifive.HiC('G1e_rep1_.hcp')
hic_2 = hifive.HiC('G1e_rep2_.hcp')

chrom = 'chr19'

'''Goal is to use this:
hifive.hic_binning.find_cis_subregion_signal(hic, chrom, binsize=10000, binbounds1=None, binbounds2=None, start1=None, stop1=None, startfend1=None, stopfend1=None, start2=None, stop2=None, startfend2=None, stopfend2=None, datatype='enrichment', skipfiltered=False, returnmapping=False, \**kwargs)'''

#First I want to store the data for loc of paired before selection ~/taylorLab/VISION/data/eRP_scores_PSU_regression/mouse_chr19_allInfo.txt
TSS_pairs = {}
all_paired_file, interaction_length = sys.argv[1:]
for line in open(all_paired_file):
    fields = line.strip('\r\n').split('\t')
    TSS = fields[1]
    if TSS not in TSS_pairs:
        TSS_pairs[TSS] = {"ccREs": [],
                          "sels": []}
    UCI = fields[8] #unique ccRE identifier chr_start_stop
    TSS_pairs[TSS]["ccREs"].append(UCI)
    sel = fields[6]
    TSS_pairs[TSS]["sels"].append(sel)


#Next I want to transfer this information to be in bin sort of space where the TSS is centered in a 2kb bin and the ccRE is centered in a 2kb bin to get the start and stop locations to compare
interaction_length = int(interaction_length)
interaction_frequencies_1 = []
interaction_frequencies_2 = []
times_selected = []
flag1 = 0
flag1s = []
flag2 = 0
flag2s = []
for key in TSS_pairs:
    start1 = int(key) - interaction_length//2
    stop1 = int(key) + interaction_length//2
    for CRE_pair, sel in zip(TSS_pairs[key]["ccREs"], TSS_pairs[key]["sels"]):
        chr, cre_start, cre_stop = CRE_pair.split('_')
        cre_mid = (int(cre_start) + int(cre_stop))//2
        start2 = cre_mid - interaction_length//2
        stop2 = cre_mid + interaction_length//2
        freq1 = hifive.hic_binning.find_cis_subregion_signal(hic_1, chrom, binsize=interaction_length, start1=start1, stop1=stop1, start2=start2, stop2=stop2, datatype='fend')
        freq2 = hifive.hic_binning.find_cis_subregion_signal(hic_2, chrom, binsize=interaction_length, start1=start1, stop1=stop1, start2=start2, stop2=stop2, datatype='fend')
        obs1 = freq1[0,0,0]
        exp1 = freq1[0,0,1]
        if exp1 == 0:
            flag1 += 1
            flag1s.append(obs1)
            exp1 = 0.001
        obs2 = freq2[0,0,0]
        exp2 = freq2[0,0,1]
        if exp2 == 0:
            flag2 += 1
            flag2s.append(obs2)
            exp2 = 0.001
        interaction1 = obs1/exp1
        interaction_frequencies_1.append(interaction1)
        interaction2 = obs2/exp2
        interaction_frequencies_2.append(interaction2)
        if sel == '0' or sel == '.':
            times_selected.append(0)
        else:
            times_selected.append(int(sel))

print(flag1, flag2)
print(np.sum(flag1s != 0), np.sum(flag2s != 0))
interaction_frequencies_1 = np.array(interaction_frequencies_1)
interaction_frequencies_2 = np.array(interaction_frequencies_2)
times_selected = np.array(times_selected)

datasets = [interaction_frequencies_1, interaction_frequencies_2]
labels = ['r1: 0-11', 'r2: 0-11']
datasets_goe = [interacion_frequencies_1,interaction_frequencies_2]
labels_goe = ['r1: 0-11', 'r2: 0-11']
for i in range(0, np.amax(times_selected)+1):
    print("i: " ,i)
    mask = np.where(times_selected == i)[0]
    datasets.append(interaction_frequencies_1[mask])
    labels.append('r1: {}'.format(i))
    datasets.append(interaction_frequencies_2[mask])
    labels.append('r2: {}'.format(i))
    mask_goe = np.where(times_selected >= i)[0]
    datasets_goe.append(interaction_frequencies_1[mask_goe])
    labels_goe.append('r1: >={}'.format(i))
    datasets_goe.append(interaction_frequencies_2[mask_goe])
    labels_goe.append('r2: >={}'.format(i))


fig, ax = plt.subplots(figsize=(40,20))
ax.set_title('Chr19 G1e Interaction Frequencies')
ax.violinplot(datasets, positions=np.arange(len(datasets)))
ax.set_ylabel("Enrichment Interaction Frequency (observed/expected)")
ax.set_xlabel("Number of times selected for inclusion")
ax.set_xticks(np.arange(len(datasets)))
ax.set_xticklabels(labels)
fig.savefig("chr19_G1e_int_freq_num.png")
plt.close(fig)

fig, ax = plt.subplots(figsize=(40,20))
ax.set_title('Chr19 G1e Interaction Frequencies')
ax.violinplot(datasets_goe, positions=np.arange(len(datasets_goe)))
ax.set_ylabel("Enrichment Interaction Frequency (observed/expected)")
ax.set_xlabel("Number of times selected for inclusion")
ax.set_xticks(np.arange(len(datasets_goe)))
ax.set_xticklabels(labels_goe)
fig.savefig("chr19_G1e_int_freq_goe.png")
plt.close(fig)
