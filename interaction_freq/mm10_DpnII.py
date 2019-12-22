#!/usr/bin/env python3

'''Adapted from Mike's Script which Snakefile calls on Comp'''
import fasta
import numpy as np
import re
import sys
genome = sys.argv[1]

RE_site = 'GATC'

genome_reader = fasta.FASTAReader(open(genome))

output_file = open('mm10_DpnII_frag.bed', 'w+')
fragments = {}
for ident, sequence in genome_reader:
    matches = np.array([match.start() for match in re.finditer(RE_site, sequence.upper())], dtype=np.int32)
    fragments[ident] = np.empty((matches.shape[0] -1, 2), dtype=np.int32)
    fragments[ident][:,0] = matches[:-1]
    fragments[ident][:,1] = matches[1:]
    for j in range(fragments[ident].shape[0]):
        output_file.write('{}\t{}\t{}\n'.format(ident, fragments[ident][j,0], fragments[ident][j,1]))

output_file.close()
