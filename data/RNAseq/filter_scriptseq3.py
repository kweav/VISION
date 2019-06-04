#!/usr/bin/env python3

import argparse as ap
import numpy as np

'''
Usage: ./filter_scriptseq3.py --TPM_file scriptseq3.v3.kw1.tab --geneLoc scriptseq3.v3.mart_export.txt
'''
parser = ap.ArgumentParser(description='filter and parse TPM files to produce a file with gene location and TPM values for each cell type')
parser.add_argument('--TPM_file',action= 'store', nargs=1, type=str, required=True, help='List of TPM files: cellType.amit.geneID_TPM.tsv')
parser.add_argument('--geneLoc', action='store', nargs=1, type=str, required=True, help='List of gene Loc files: cellType.amit.mart_export.txt')
args = parser.parse_args()
TPM_file = args.TPM_file[0]
geneLoc_file = args.geneLoc[0]

geneID_to_loc = {}
for i, line in enumerate(open(geneLoc_file)):
    if i == 0:
        pass
    else:
        fields = line.strip('\r\n').split('\t')
        geneID = fields[0]
        chr = fields[1]
        start = fields[2]
        end = fields[3]
        if geneID not in geneID_to_loc:
            geneID_to_loc[geneID] = (chr, start, end)

geneID_cellType_TPM = {}
field_to_cellType = {2: 'LSK',
                     3: 'LSK',
                     4: 'CMP',
                     5: 'CMP',
                     6: 'GMP',
                     7: 'GMP',
                     8: 'MEP',
                     9: 'MEP',
                     10: 'CFU-E',
                     11: 'CFU-E',
                     12: 'ERYad',
                     13: 'ERYad',
                     14: 'CFU-Mk',
                     15: 'CFU-Mk',
                     16: 'iMK',
                     17: 'iMK',
                     18: 'Mono',
                     19: 'Mono',
                     20: 'Neutrophil',
                     21: 'Neutrophil',
                     22: 'G1E',
                     23: 'G1E',
                     24: 'ER4',
                     25: 'ER4'}
cellTypes = ['LSK', 'CMP', 'GMP', 'MEP', 'CFU-E', 'ERYad', 'CFU-Mk', 'iMK', 'Mono', 'Neutrophil', 'G1E', 'ER4']
for line in open(TPM_file):
    fields=line.strip('\r\n').split('\t')
    geneID = fields[0].split(".")[0]
    if geneID not in geneID_cellType_TPM:
        geneID_cellType_TPM[geneID] = {}
    for i in range(2,26):
        tpm = float(fields[i])
        cellType = field_to_cellType[i]
        if cellType not in geneID_cellType_TPM[geneID]:
            geneID_cellType_TPM[geneID][cellType] = []
        geneID_cellType_TPM[geneID][cellType].append(tpm)

file=open("scriptseq3.v3.kw2.tab", "w+")
for key in geneID_cellType_TPM:
    if key not in geneID_to_loc:
        print(key)
    else:
        chr, start, stop = geneID_to_loc[key]
        if chr == 'MT' or 'PATCH' in chr:
            pass
        else:
            file.write(chr + '\t' + start + '\t' + stop)
            for cellType in cellTypes:
                tpmValue = np.average(geneID_cellType_TPM[key][cellType])
                file.write('\t' + str(tpmValue))
            file.write('\n')
