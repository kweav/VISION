#!/usr/bin/env python3

import argparse as ap
import numpy as np

parser = ap.ArgumentParser(description='retrieve genomic locations for geneIDs from RNAseq data using Gencode annotation')
parser.add_argument('--gencode', action='store', nargs=1, type=str, required=True, help='Gencode annotation')
parser.add_argument('--RNAseq', action='store', nargs = '+', type=str, required=True, help='RNAseq files with geneIDs')
args=parser.parse_args()
gencode_file = args.gencode[0]
RNAseq_files = args.RNAseq

geneID_to_loc = {}
for line in open(gencode_file):
    if line.startswith('#'):
        pass
    else:
        fields=line.strip('\r\n').split('\t')
        if fields[2] == 'gene':
            chr = fields[0]
            start = fields[3]
            stop = fields[4]
            info = fields[8]
            geneID = info.split(";")[0].split(" ")[1].replace('"','')
            if geneID not in geneID_to_loc:
                geneID_to_loc[geneID] = (chr, start, stop)

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
cellTypes2 = ['B', 'NK', 'T_CD4', 'T_CD8']

geneID_cellType_TPM = {}
geneID_cellType_TPM1 = {}
for RNAseq_file in RNAseq_files:
    if RNAseq_file == 'scriptseq3.v3.kw1.tab':
        for line in open(RNAseq_file):
            fields=line.strip('\r\n').split('\t')
            geneID = fields[0]
            if geneID not in geneID_cellType_TPM1:
                geneID_cellType_TPM1[geneID] = {}
                for i in range(2,26):
                    tpm = float(fields[i])
                    cellType = field_to_cellType[i]
                    if cellType not in geneID_cellType_TPM1[geneID]:
                        geneID_cellType_TPM1[geneID][cellType] = []
                    geneID_cellType_TPM1[geneID][cellType].append(tpm)

        file=open("scriptseq3.v3.kw2.tab", "w+")
        file2=open('noCorrespondingLoc_scriptseq.txt', "w+")
        for key in geneID_cellType_TPM1:
            if key not in geneID_to_loc:
                file2.write(key + '\n')
            else:
                chr, start, stop = geneID_to_loc[key]
                if 'M' in chr or 'PATCH' in chr:
                    pass
                else:
                    newFinalField = ''
                    for cellType in cellTypes:
                        tpmValue = np.average(geneID_cellType_TPM1[key][cellType])
                        newFinalField += cellType
                        newFinalField += '='
                        newFinalField += "{0:.2f}".format(tpmValue)
                        newFinalField += ";"
                    file.write(chr + '\t' + start + '\t' + stop + '\t' + newFinalField + '\n')
        file.close()
        file2.close()
    else:
        cellType_TPM = RNAseq_file.split(".")[0].replace("amitData/", "")
        for line in open(RNAseq_file):
            fields=line.strip('\r\n').split('\t')
            geneID = fields[0]
            tpm = fields[1]
            if geneID not in geneID_cellType_TPM:
                geneID_cellType_TPM[geneID] = {}
            if cellType_TPM not in geneID_cellType_TPM[geneID]:
                geneID_cellType_TPM[geneID][cellType_TPM] = 0
            geneID_cellType_TPM[geneID][cellType_TPM] = tpm

file=open("amit.cellTypes_withLoc.tab", "w+")
file2=open("noCorrespondingLoc_amit.txt", "w+")
for key in geneID_cellType_TPM:
    if key not in geneID_to_loc:
        file2.write(key + '\n')
    else:
        chr, start, stop = geneID_to_loc[key]
        if 'M' in chr or 'PATCH' in chr:
            pass
        else:
            newFinalField = ''
            for cellType in cellTypes2:
                tpmValue = geneID_cellType_TPM[key][cellType]
                newFinalField += cellType
                newFinalField += '='
                newFinalField += "{0:.2f}".format(float(tpmValue))
                newFinalField += ';'
            file.write(chr + '\t' + start + '\t' + stop + '\t' + newFinalField + '\n')
file.close()
file2.close()
