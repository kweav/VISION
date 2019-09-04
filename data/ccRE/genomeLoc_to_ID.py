#!/usr/bin/env python3

#import sys
#Usage when checking accuracy of encoding/decoding: ./genomeLoc_to_ID.py ccRE_locs.txt > IDs.txt
#Usage when producing a file with encoded locations: ./genomeLoc_to_ID.py --input_file genes_ccRE_1Mb_diff_1kb.bed --output_file genes_ccRE_ID_1Mb_diff_1kb.bed
#                                                    ./genomeLoc_to_ID.py --input_file genes_ccRE_1kb_window.bed    --output_file genes_ccRE_ID_1kb_window.bed
#                                                    ./genomeLoc_to_ID.py --input_file genes_ccRE_200bp_window.bed --output_file genes_ccRE_ID_200bp_window.bed
#                                                    ./genomeLoc_to_ID.py --input_file VISIONmusHem_ccREs_filterkw.bed --ccRE_chr_field 0 --ccRE_start_field 1 --ccRE_end_field 2 --output_file VISIONmusHem_ccREs_filterkw2.bed, then cut -f 6,7,8 

import argparse as ap

'''produces 19 character unique IDs where
    -   the first character points to chr
    -   the rest of the ID can be split on 'h.' to get the start and ending genome locations in a pseudo hexadecimal form; used 'h.' to separate start and end since h alone could be true hexadecimal
        -   First, replace all filler '.' with ''; these fillers were used so that all IDs were the same length
        -   Second, add '0x' to the front to produce the true hexadecimal
        -   Third, convert hexadecimal back to int
        -   Decoding script is ID_to_genomeLoc.py'''

#file = sys.argv[1]

parser = ap.ArgumentParser(description='encode ccREs in windows to a unique identifier')
parser.add_argument('--input_file', action='store', nargs=1, type=str, required = True, help="file with genes/expression information and the ccREs within a given window (product of bedtools window)")
parser.add_argument('--ccRE_chr_field', action='store', nargs=1, type=int, required=False, default =[5], help='field with 0 indexing where the ccRE chromosome is located')
parser.add_argument('--ccRE_start_field', action='store', nargs=1, type=int, required=False, default=[6], help='field with 0 indexing where the ccRE start is located')
parser.add_argument('--ccRE_end_field', action='store', nargs=1, type=int, required=False, default=[7], help='field with 0 indexing where the ccRE end is located')
parser.add_argument('--output_file', action='store', nargs=1, type=str, required=True, help='desired output file name where the ccRE genomic location fields are collapsed to a single field with ID')
args=parser.parse_args()
file = args.input_file[0]
l = args.ccRE_chr_field[0]
m = args.ccRE_start_field[0]
n = args.ccRE_end_field[0]
out_file = open(args.output_file[0], 'w+')

chrEncode = {'chr1':'A',
            'chr2':'B',
            'chr3':'C',
            'chr4':'D',
            'chr5':'E',
            'chr6':'F',
            'chr7':'G',
            'chr8':'H',
            'chr9':'I',
            'chr10':'J',
            'chr11':'K',
            'chr12':'L',
            'chr13':'M',
            'chr14':'N',
            'chr15':'O',
            'chr16':'P',
            'chr17':'Q',
            'chr18':'R',
            'chr19':'S',
            'chrX':'T',
            'chrY':'U'}

#lensOID = []
for line in open(file):
    ID = ''
    fields=line.strip('\r\n').split('\t')
    tot_fields = len(fields)
    chr = fields[l]
    start = int(fields[m])
    end = int(fields[n])
    ID += chrEncode[chr]
    h_start = hex(start)[2:]
    if len(h_start) < 7:
        for i in range(7 - len(h_start)):
            h_start = "."+h_start
    h_start = "h."+h_start
    ID+= h_start
    h_end = hex(end)[2:]
    if len(h_end) < 7:
        for i in range(7 - len(h_end)):
            h_end = "."+h_end
    h_end = "h."+h_end
    ID += h_end
    #lensOID.append(len(ID))
    lineToWrite = ''
    for j in range(tot_fields-l):
        lineToWrite += fields[j]
        lineToWrite += '\t'
    lineToWrite += ID
    lineToWrite += '\t'
    for k in range(n+1,tot_fields):
        lineToWrite += fields[k]
        if k != tot_fields - 1:
            lineToWrite += '\t'
        else:
            lineToWrite += '\n'
    out_file.write(lineToWrite)
#print(set(lensOID))
