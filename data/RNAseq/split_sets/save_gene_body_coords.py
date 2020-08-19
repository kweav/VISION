#!/usr/bin/env python3

import sys
import numpy as np

#usage:
'''
./save_gene_body_coords.py rnaTPM_train_pc.txt 19012 trainGeneCoords_pc.npz
./save_gene_body_coords.py rnaTPM_test_pc.txt 1967 testGeneCoords_pc.npz
./save_gene_body_coords.py rnaTPM_ref_pc.txt 857 refGeneCoords_pc.npz
'''

input_file = sys.argv[1]
tssN = int(sys.argv[2])

to_save = np.zeros((tssN, 3), dtype=np.object)

for i, line in enumerate(open(input_file)):
    fields = line.strip('\r\n').split()
    if i != 0:
        to_save[i-1,0] = fields[0] #set chr
        to_save[i-1,1] = int(fields[1]) #set "start", but really just the minimum gene body edge value
        to_save[i-1,2] = int(fields[2]) #set "end", but really just the maximum gene body edge value

output_file = sys.argv[3]
f = open(output_file, 'wb')
np.savez(f, geneBodyCoords = to_save)
f.close()
