#!/usr/bin/env python3

import argparse as ap
import numpy as np

'''Usage:./geneID_loc_exp.py --gencode ~/mm10_genome/gencode.vM4.annotation.gtf --RNAseq [scriptseq3.v3.kw1.tab] or [amitData/*.amit.geneID_TPM.tab]'''

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
field_to_cellType = {2: 'Lsk',
                     3: 'Lsk',
                     4: 'Cmp',
                     5: 'Cmp',
                     6: 'Gmp',
                     7: 'Gmp',
                     8: 'Mep',
                     9: 'Mep',
                     10: 'Cfue',
                     11: 'Cfue',
                     12: 'Eryad',
                     13: 'Eryad',
                     14: 'Cfum',
                     15: 'Cfum',
                     16: 'Imk',
                     17: 'Imk',
                     18: 'Mon',
                     19: 'Mon',
                     20: 'Neu',
                     21: 'Neu',
                     22: 'G1e',
                     23: 'G1e',
                     24: 'Er4',
                     25: 'Er4'}
cellTypes = ['Lsk', 'Cmp', 'Gmp', 'Mep', 'Cfue', 'Eryad', 'Cfum', 'Imk', 'Mon', 'Neu', 'G1e', 'Er4']
cellTypes2 = ['B', 'Nk', 'Tcd4', 'Tcd8']

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

        file=open("scriptseq3.v3.kw2.IDlocexp.bed", "w+")
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
                    file.write(chr + '\t' + start + '\t' + stop + '\t' + key + '\t' + newFinalField + '\n')
        file.close()
        file2.close()
    else:
        cellType_TPM = RNAseq_file.split(".")[0].replace("amitData/", "")
        print(cellType_TPM)
        for line in open(RNAseq_file):
            fields=line.strip('\r\n').split('\t')
            geneID = fields[0]
            tpm = fields[1]
            if geneID not in geneID_cellType_TPM:
                geneID_cellType_TPM[geneID] = {}
            if cellType_TPM not in geneID_cellType_TPM[geneID]:
                geneID_cellType_TPM[geneID][cellType_TPM] = 0
            geneID_cellType_TPM[geneID][cellType_TPM] = tpm

if RNAseq_file != 'scriptseq3.v3.kw1.tab':
    file=open("amit.cellTypes_IDlocexp.bed", "w+")
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
                file.write(chr + '\t' + start + '\t' + stop  + '\t' + key + '\t' +  newFinalField + '\n')
    file.close()
    file2.close()
