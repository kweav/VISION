#!/usr/bin/env python3

'''Usage: ./prepare_expression.py rnaTPM.txt 41814 TSS_expression.npz
          ./prepare_expression.py rnaTPM_pc_only.txt 21782 TSS_pc_expression.npz 1'''

import sys
import numpy as np
import os

input_file = sys.argv[1]
tssN = int(sys.argv[2])
if len(sys.argv) > 4:
    threads = sys.argv[4]
else:
    threads = 1

os.environ["OMP_NUM_THREADS"] = threads
os.environ["OPENBLAS_NUM_THREADS"] = threads
os.environ["MKL_NUM_THREADS"] = threads
os.environ["VECLIB_MAXIMUM_THREADS"] = threads
os.environ["NUMEXPR_NUM_THREADS"] = threads


header_to_index = {"CFU_E_ad": "Cfue",
                   "CFUMK": "Cfum",
                   "CMP": "Cmp",
                   "ERY_fl": "Eryfl",
                   "GMP": "Gmp",
                   "MK_imm_ad": "Imk",
                   "LSK_BM": "Lsk",
                   "MEP": "Mep",
                   "MONO_BM": "Mon",
                   "NEU": "Neu",
                   "ER4": "Er4",
                   "G1E": "G1e"}
exp = np.zeros((tssN, len(header_to_index.keys())+2), dtype=np.object)
for i, line in enumerate(open(input_file)):
    fields = line.strip('\r\n').split()
    if i == 0:
        header_cellTypes = fields[4:]
    else:
        exp[i-1, 0:2] = fields[0:2] #set chr and TSS
        exp[i-1, 2:] = fields[4:] #set expression

cellIndex = np.empty(len(header_cellTypes), dtype=np.object)
for j, ct in enumerate(header_cellTypes):
    cellIndex[j] = header_to_index[ct]

output_file = sys.argv[3]
f = open(output_file, 'wb')
np.savez(f, exp = exp[:, 2:].astype(np.float32), TSS = exp[:, 0:2], cellIndex = cellIndex)
f.close()
