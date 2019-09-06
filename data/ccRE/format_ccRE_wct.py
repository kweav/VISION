#!/usr/bin/env python3

import sys

"""
Usage: ./format_ccRE_wct.py VISIONmusHem_ccREs.bed > VISIONmusHem_ccREs_wct_filterkw.bed
"""

cellTypesIS = ["Lsk", "Hpc7", "Cmp", "Mep", "G1e", "Er4", "Cfue", "Eryad", "Eryfl", "Cfum", "Imk", "Gmp", "Mon", "Neu", "Nk", "B", "Tcd4", "Tcd8"]
cellTypesSA = ["Lsk", "Hpc7", "Cmp", "Mep", "G1e", "Er4", "Cfue", "Eryad", "Eryfl", "Cfum", "Imk", "Mk", "Gmp", "Mon", "Neu", "Clp", "Nk", "B", "Tcd4", "Tcd8"]

for line in open(sys.argv[1]):
    fields=line.strip("\r\n").split("\t")
    newLine1 = ""
    minifield = fields[3].split("_")
    for j, cellTypeIS in zip(range(len(minifield)), cellTypesIS):
        addition = cellTypeIS + '=' + minifield[j] + ';'
        newLine1 += addition
    newLine2=""
    for i,cellTypeSA in zip(range(4,4+len(fields[4:])),cellTypesSA):
        addition = cellTypeSA + '=' + fields[i] + ';'
        newLine2 += addition
    print(fields[0], fields[1], fields[2], newLine1, newLine2, sep='\t')
