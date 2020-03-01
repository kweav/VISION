#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd

chrom = sys.argv[1]
interactions_to_query = sys.argv[2]
header_line = sys.argv[3]
sparse_num = sys.argv[4] #2048?

cellTypes = ['B', 'CD4', 'CD8', 'G1E', 'HPC7', 'LSK', 'Neu']

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

partition_file_base = '/project/vision/HiC_Modeling/HP_Annotations/mm10_{}_{}.npy' #first {} for restriction enzyme #second {} for chrom
sparse_base = '/project/vision/HiC_Modeling/Binned_HiC/{}/{}_{}_{}.npy' #first and second {} for CT. third {} for sparse number and 4th {} for chrom
restriction_enzyme_file = '/project/vision/HiC_Modeling/datasets.txt'
npz_base = '/project/vision/Target_Genes/analysis_PSU_method/int_freq_{}_{}.npz' #first {} for chrom second {} for partition of that chrom
bias_file_base = '/project/vision/HiC_Modeling/Binned_HiC/{}/Bias/{}_{}_{}_bias.npy'  #first and second {} for CT. third {} for sparse number and 4th {} for chrom

cellType_to_RE = {}
for line in open(restriction_enzyme_file):
    fields = line.strip('\r\n').split()
    CT = fields[1]
    RE = fields[2]
    if CT not in cellType_to_RE:
        cellType_to_RE[CT] = RE

RE_locs_to_bins = {}
for RE in np.unique(list(cellType_to_RE.values())):
    if RE != "CviQI-CbiAII-BfaI":
        RE_locs_to_bins[RE] = np.zeros(chrSizes[chrom])
        bin_annotations = np.load(partition_file_base.format(RE, chrom))
        for i in range(bin_annotations.shape[0]):
            start, stop = bin_annotations[i, :]
            RE_locs_to_bins[RE][start:stop] = i

'''read in interaction to query as a pandas dataframe'''
def read_in(interactions_to_query, header_line):
    df_int = pd.read_csv(interactions_to_query, sep='\t', header=None)
    df_int.columns = open(header_line).readline().strip('\r\n').split('\t')
    df_int_sub = df_int.iloc[:, 0:9]
    return(df_int_sub)

df_int_sub = read_in(interactions_to_query, header_line)
tss_unique = np.unique(df_int_sub['tss'])
tssN = tss_unique.shape[0]
tss_to_index = {tss:i for i, tss in enumerate(tss_unique)}
tss_indices = list(map(lambda x: tss_to_index[x], np.array(df_int_sub['tss'])))
cre_unique = np.unique(df_int_sub['pk'])
creN = cre_unique.shape[0]
cre_to_index = {cre:i for i, cre in enumerate(cre_unique)}
cre_indices = list(map(lambda x: cre_to_index[x], np.array(df_int_sub['pk'])))

num_selected = np.full((tssN, creN), np.nan)
df_int_sub['selected'].replace(to_replace = '.', value=0, inplace=True)
num_selected[tss_indices, cre_indices] = np.array(df_int_sub['selected']).astype(np.int32) #USE THESE ARRAYS TO PLOT DISTRIBUTION OF NUM PAIRED FOR EACH GENE!

cre_mid = np.array([(x+y)/2 for x,y in zip(df_int_sub['ccRE_start'], df_int_sub['ccRE_end'])]).astype(np.int32)
int_freq = np.full((tssN, creN, len(cellTypes)), np.nan)
for j, CT in enumerate(cellTypes):
    RE = cellType_to_RE[CT]
    biases = np.load(bias_file_base.format(CT, CT, sparse_num, chrom))
    cre_mid_bins = RE_locs_to_bins[RE][cre_mid].astype(np.int32)
    tss_bins = RE_locs_to_bins[RE][list(df_int_sub['tss'])].astype(np.int32)
    concat = np.concatenate((cre_mid_bins.reshape((-1, 1)), tss_bins.reshape((-1,1))), axis=1)

    bin1 = np.amin(concat, axis=1)
    bin1_arg = np.argmin(concat, axis=1)
    bin2 = np.amax(concat, axis=1)
    bin2_arg = np.argmax(concat, axis=1)

    sparse_binned = np.load(sparse_base.format(CT, CT, sparse_num, chrom))


    bin1_bool = np.in1d(bin1, sparse_binned['bin1'])
    int_freq[np.array(tss_indices)[~bin1_bool], np.array(cre_indices)[~bin1_bool],j] = 0

    vals, first_indexes, counts = np.unique(sparse_binned['bin1'], return_index=True, return_counts=True)

    corresponding_bin2_sparse = np.full((vals.shape[0], np.max(counts)), np.nan)
    val_to_corresp_index = {}
    for k, (val, index, count) in enumerate(zip(vals, first_indexes, counts)):
        if k == first_indexes.shape[0]-1:
            corresponding_bin2_sparse[k, 0:count] = sparse_binned['bin2'][index:]
        else:
            corresponding_bin2_sparse[k, 0:count] = sparse_binned['bin2'][index: first_indexes[k+1]]
        val_to_corresp_index[val] = (k, index)

    nonzero_bin1 = bin1[bin1_bool]
    nonzero_bin2 = bin2[bin1_bool]
    nonzero_tss_indices = np.array(tss_indices)[bin1_bool]
    nonzero_cre_indices = np.array(cre_indices)[bin1_bool]

    for l, (nzb1, nzb2) in enumerate(zip(nonzero_bin1, nonzero_bin2)):
        k, first_index = val_to_corresp_index[nzb1]
        #get the score
        if np.in1d(nzb2, corresponding_bin2_sparse[k,:]):
            location = np.where(corresponding_bin2_sparse[k,:] == nzb2)[0]
            score_val = sparse_binned['score'][location+first_index]
            #adjust the score
            bias_adj = biases[nzb1]*biases[nzb2]
            if bias_adj == 0:
                adj_score = 0
            else:
                adj_score = score_val/bias_adj
        else:
            adj_score = 0
        #store the score
        int_freq[np.array(nonzero_tss_indices[l]), np.array(nonzero_cre_indices[l]), j] = adj_score

outfile = 'original_paired_{}.npz'.format(chrom)
f = open(outfile, 'wb')
np.savez(f, int_freq = int_freq, num_selected = num_selected)
f.close()
