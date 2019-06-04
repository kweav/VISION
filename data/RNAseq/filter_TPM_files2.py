#!/usr/bin/env python3

import argparse as ap

'''
Usage: ./filter_TPM_files2.py --TPM_files amit.geneID_TPM.tsvFiles --geneLoc amit.mart_export.txtFiles
'''
parser = ap.ArgumentParser(description='filter and parse TPM files to produce a file with gene location and TPM values for each cell type')
parser.add_argument('--TPM_files',action= 'store', nargs='+', type=str, required=True, help='List of TPM files: cellType.amit.geneID_TPM.tsv')
parser.add_argument('--geneLoc', action='store', nargs='+', type=str, required=True, help='List of gene Loc files: cellType.amit.mart_export.txt')
args = parser.parse_args()
TPM_file_list = args.TPM_files
geneLoc_file_list = args.geneLoc

geneID_to_loc = {}
for geneLoc_file in geneLoc_file_list:
    cellType_geneLoc = geneLoc_file.split(".")[0].replace("amitData/", "")
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
for TPM_file in TPM_file_list:
    cellType_TPM = TPM_file.split(".")[0].replace("amitData/", "")
    for line in open(TPM_file):
        fields=line.strip('\r\n').split('\t')
        geneID = fields[0]
        tpm = fields[1]
        if geneID not in geneID_cellType_TPM:
            geneID_cellType_TPM[geneID] = {}
        if cellType_TPM not in geneID_cellType_TPM[geneID]:
            geneID_cellType_TPM[geneID][cellType_TPM] = 0
        geneID_cellType_TPM[geneID][cellType_TPM] = tpm

file=open("amit.cellTypes_withLoc.tsv", "w+")
for key in geneID_cellType_TPM:
    if key not in geneID_to_loc:
        print(key)
    else:
        chr, start, stop = geneID_to_loc[key]
        if chr == 'MT' or 'PATCH' in chr:
            pass
        else:
            file.write(chr + '\t' + start + '\t' + stop)
            for cellType in ['B', 'NK', 'T_CD4', 'T_CD8']:
                tpmValue = geneID_cellType_TPM[key][cellType]
                file.write('\t' + tpmValue )
            file.write('\n')
