#!/usr/bin/env python3

import sys

"""
Usage: ./format_ccRE.py VISIONmusHem_ccREs.bed > VISIONmusHem_ccREs_filterkw.bed
"""
for line in open(sys.argv[1]):
    fields=line.strip("\r\n").split("\t")
    newLine=""
    for i in range(4,4+len(fields[4:])):
        newLine += (fields[i]+";")
    print(fields[0], "\t", fields[1], "\t", fields[2], "\t", fields[3], '\t', newLine, sep='')
