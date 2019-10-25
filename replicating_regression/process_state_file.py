#!/usr/bin/env python3

import numpy as np

statepref = '/home/kweave23/VISION_regression/their_stuff/pknorm_2_16lim_ref1mo_0424_lesshet'
exp_file = '/home/kweave23/VISION_regression/their_stuff/rnaTPM.txt'
chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
                'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
                'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
def load_RNA(exp_file): #rnaTPM.txt
    with open(exp_file) as f:
        rna_names = f.readline().split()[4:]
    return (rna_names)

def load_state(statepref, chr, rna_names): #pknorm_2_16lim_ref1mo_0424_lesshet.state
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
    t = np.where(scr == chr)[0] #extract only relevant chromosomes genes/states;
    pos = pos[t] #extract only relevant chromosome state bins
    state = state[t,:] #extract only relevant chromosomes states
    return (state, pos)

rna_names = load_RNA(exp_file)
to_save = {}
f = open('state_and_pos_by_chr.npz', 'wb')
for chr in chromosomes:
    state, pos = load_state(statepref, chr, rna_names)
    to_save["{}_state".format(chr)] = state
    to_save["{}_pos".format(chr)] = pos
np.savez(f, **to_save)
f.close()
