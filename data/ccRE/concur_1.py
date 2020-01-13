#!/usr/bin/env python3

'''Usage: ./concur_1.py VISIONmusHem_ccREs_wct_filterkw_prop.bed'''

import sys

for line in open(sys.argv[1]):
    fields=line.strip('\r\n').split('\t')
    subfields = fields[4].split(';')[:-1]
    for sub in subfields:
        cellType, props = sub.split('=')
        splits = props.split(')')[:-1]
        if len(splits) == 1:
            state, prop = splits[0].split('_')
            if float(prop) != 1.0:
                print("Error A: ", prop)
        if len(splits) > 1:
            full_prop = 0
            states = []
            for partial in splits:
                state, prop = partial.split('_')
                states.append(state.replace('(', ''))
                full_prop += float(prop)
            state_set = set(states)
            if len(state_set) != len(splits):
                print("Error B")
            if round(full_prop) != 1.0:
                print('Error C: ', full_prop)
