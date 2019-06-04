#!/usr/bin/env python3

import argparse as ap

'''
Usage: ./filter_TPM_files.py --TPM_files amit.tsvFiles
'''
parser = ap.ArgumentParser(description='filter and parse TPM files to produce a file with geneID and TPM value and a second file of just geneID for biomart')
parser.add_argument('--TPM_files',action= 'store', nargs='+', type=str, required=True, help='List of TPM files')
args = parser.parse_args()
TPM_file_list = args.TPM_files

for TPM_file in TPM_file_list:
    cellType = TPM_file.split(".")[0]

    f1 = open("{}.amit.geneID_TPM.tsv".format(cellType), "w+")
    f2 = open("{}.amit.geneID.txt".format(cellType), "w+")
    for line in open(TPM_file):
        if line.startswith("E"):
            fields = line.strip('\r\n').split('\t')
            geneID, suffix = fields[0].split(".")
            TPM = fields[5]
            f1.write(geneID + '\t' + TPM + '\n')
            f2.write(geneID + '\n')

    f1.close()
    f2.close()
