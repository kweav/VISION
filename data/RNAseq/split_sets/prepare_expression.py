#!/usr/bin/env python3

'''Usage: ./prepare_expression.py input_file length output_npz'''

import sys
import numpy as np

input_file = sys.argv[1]
tssN = int(sys.argv[2])

for i, line in enumerate(open(input_file)):
    fields = line.strip('\r\n').split()
    if i == 0:
        header_cellTypes = fields[6:]
        exp = np.zeros((tssN, len(header_cellTypes)+2), dtype=np.object)
    else:
        #Find TSS
        if fields[5] == '+':
            TSS = fields[1]
        elif fields[5] == '-':
            TSS = fields[2]
        #Store info
        exp[i-1, 0] = fields[0] #set chr
        exp[i-1, 1] = TSS #set TSS
        exp[i-1, 2:] = fields[6:] #set expression

cellIndex = np.empty(len(header_cellTypes), dtype=np.object)
for j, ct in enumerate(header_cellTypes):
    cellIndex[j] = ct

output_file = sys.argv[3]
f = open(output_file, 'wb')
np.savez(f, exp = exp[:, 2:].astype(np.float32), TSS = exp[:, 0:2], cellIndex = cellIndex)
f.close()
