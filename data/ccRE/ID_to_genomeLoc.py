#!/usr/bin/env python3

import sys

'''Usage: ./ID_to_genomeLoc.py IDs.txt > decode_locs.bed
decodes 19 character unique IDs where
    -   the first character points to chr
    -   the rest of the ID can be split on 'h.' to get the start and ending genome locations in a pseudo hexadecimal form; used 'h.' to separate start and end since h alone could be true hexadecimal
        -   First, replace all filler '.' with ''; these fillers were used so that all IDs were the same length
        -   Second, add '0x' to the front to produce the true hexadecimal
        -   Third, convert hexadecimal back to int
        -   ID producing script is genomeLoc_to_ID.py'''

chrDecode = {'A':'chr1',
            'B':'chr2',
            'C':'chr3',
            'D':'chr4',
            'E':'chr5',
            'F':'chr6',
            'G':'chr7',
            'H':'chr8',
            'I':'chr9',
            'J':'chr10',
            'K':'chr11',
            'L':'chr12',
            'M':'chr13',
            'N':'chr14',
            'O':'chr15',
            'P':'chr16',
            'Q':'chr17',
            'R':'chr18',
            'S':'chr19',
            'T':'chrX',
            'U':'chrY'}

def decode_string(chrDecode, ID):
    chrD = chrDecode[ID[0]]
    filler,start,end = ID[1:].split('h.')
    start_hex = '0x'+ start.replace(".",'')
    end_hex = '0x'+ end.replace('.','')
    return(chrD, start_hex, end_hex)

file = sys.argv[1]

for line in open(file):
    lineToPrint = ''
    ID=line.strip('\r\n')
    chrD, start_hex, end_hex = decode_string(chrDecode, ID)
    lineToPrint += chrD
    lineToPrint += "\t"
    lineToPrint += str(int(start_hex,0))
    lineToPrint += '\t'
    lineToPrint += str(int(end_hex,0))
    print(lineToPrint)
