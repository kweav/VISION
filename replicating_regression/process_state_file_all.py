#!/usr/bin/env python3

import numpy as np
import pickle

statepref = '/home/kweave23/VISION_regression/their_stuff/pknorm_2_16lim_ref1mo_0424_lesshet'
exp_file = '/home/kweave23/VISION_regression/their_stuff/rnaTPM.txt'

def load_RNA(exp_file): #rnaTPM.txt
    with open(exp_file) as f:
        rna_names = f.readline().split()[4:]
    return (rna_names)

def load_state(statepref, rna_names): #pknorm_2_16lim_ref1mo_0424_lesshet.state
    file_to_load = "%s.state" %statepref
    state = []
    with open(file_to_load) as f:
        header = f.readline().split()[4:]
        for line in f:
            state.append(line.strip('\r\n').split())

    state = np.array(state, dtype=np.object)

    #state data is a file with columns for binID, chr, posStart, posEnd, celltype1, ..., cellTypeN, posClass
    scr = state[:,1] #state chromosome names
    pos = state[:,2].astype(np.int32)//200 #state bins (coordinate --> 200bp bin space)
    state = state[:,4:].astype(np.int32) #remove state position data, so just states left
    # valid = [x for x in header if x in rna_names]
    # state = state[:,np.where(valid)[0]] #select only cell types RNA data; in essence subselecting and reordering the state cell types to match the rna cell types
    valid = [header.index(x) for x in rna_names if x in header]
    state = state[:,valid] #select only cell types RNA data; in essence subselecting and reordering the state cell types to match the rna cell types

    '''subselect chromosomes'''
    chroms = np.unique(scr)
    state_all = {}
    pos_all = {}
    for chrom in chroms:
        t = np.where(scr == chrom)[0]
        pos_all[chrom] = pos[t]
        state_all[chrom] = state[t,:]

    return ( state_all, pos_all )

rna_names = load_RNA(exp_file)

with open('state_all_and_pos_all_by_chr.pickle', 'wb') as f:
    state_all, pos_all = load_state(statepref, rna_names)
    to_pickle = {'state_all': state_all,
                 'pos_all': pos_all}
    pickle.dump(to_pickle, f, protocol=pickle.HIGHEST_PROTOCOL)
