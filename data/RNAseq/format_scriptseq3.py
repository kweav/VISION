#!/usr/bin/env python3

import argparse as ap

parser = ap.ArgumentParser(description = 'format scriptseq3 file and format Amit file')
parser.add_argument('--scriptseq', action='store', nargs=1, required=True, help='scriptseq3 file with 12 cell types')
parser.add_argument('--amit', action='store', nargs=1, required=True, help='Amit data with 4 cell types')
args = parser.parse_args()
scriptseq_file = args.scriptseq[0]
amit_file = args.amit[0]

field_to_cellType = {3: 'Lsk',
                     4: 'Cmp',
                     5: 'Gmp',
                     6: 'Mep',
                     7: 'Cfue',
                     8: 'Eryad',
                     9: 'Cfum',
                     10: 'Imk',
                     11: 'Mon',
                     12: 'Neu',
                     13: 'G1e',
                     14: 'Er4'}

field_to_cellType2 = {3: 'B',
                      4: 'NK',
                      5: 'T_CD4',
                      6: 'T_CD8'}

newFile = open("scriptseq3.v3.kw3.tab", "w+")
for line in open(scriptseq_file):
    fields = line.strip('\r\n').split('\t')
    chr = fields[0]
    newChr = "chr" + chr
    start = fields[1]
    stop = fields[2]
    newFinalField = ''
    for i in range(3, 15):
        newFinalField += field_to_cellType[i]
        newFinalField += '='
        newFinalField += "{0:.2f}".format(float(fields[i]))
        newFinalField += ";"
    newFile.write(newChr + '\t' + start + '\t' + stop + '\t' + newFinalField + '\n')

newFile2 = open("amit.cellTypes_withLoc.v2.tab", "w+")
for line in open(amit_file):
    fields=line.strip('\r\n').split('\t')
    chr=fields[0]
    newChr = "chr" + chr
    start = fields[1]
    stop = fields[2]
    newFinalField = ''
    for i in range(3,7):
        newFinalField += field_to_cellType2[i]
        newFinalField += '='
        newFinalField += "{0:.2f}".format(float(fields[i]))
        newFinalField += ";"
    newFile2.write(newChr + '\t' + start + '\t' + stop + '\t' + newFinalField + '\n')
