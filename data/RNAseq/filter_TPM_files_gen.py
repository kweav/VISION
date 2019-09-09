#!/usr/bin/env python3

import argparse as ap

'''
Usage: ./filter_TPM_files.py --TPM_file scriptseq3.v3.tab
'''
parser = ap.ArgumentParser(description='filter and parse TPM files to produce a file with geneID and TPM value')
parser.add_argument('--TPM_file',action= 'store', nargs=1, type=str, required=True, help='List of TPM files')
args = parser.parse_args()
TPM_file = args.TPM_file[0]

f1 = open("scriptseq3.v3.kw1.test.tab", "w+")

for line in open(TPM_file):
    if line.startswith("E"):
        f1.write(line.strip('\r\n') + '\n')

f1.close()
