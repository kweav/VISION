#!/usr/bin/env python3

import sys
#import subprocess
import numpy as np

to_subset_file = sys.argv[1]
subset_by_file = sys.argv[2]
to_save = sys.argv[3]

protein_coding_geneIDs = []
for line in open(subset_by_file):
    fields=line.strip('\r\n').split('\t')
    subfields = fields[8].split(';')
    gene_id = subfields[0].split()[1].replace('"', '')
    gene_type = subfields[2].split()[1].replace('"', '')
    if gene_type == 'protein_coding':
        protein_coding_geneIDs.append(gene_id)
print(len(protein_coding_geneIDs))

lines_to_print = []
for line in open(to_subset_file):
    fields=line.strip('\r\n').split('\t')
    gene_id = fields[2].split(':')[-1]
    if gene_id in protein_coding_geneIDs:
        lines_to_print.append(line.strip('\r\n'))

print(len(lines_to_print))

lines_to_print = np.array(lines_to_print, dtype=np.object).reshape(-1, )
np.savetxt(to_save, lines_to_print, fmt='%s')
# #joined_lines = ','.join(lines_to_sed)
# #subprocess.call("sed -n -e{{}}'p' {} > pc_{}".format(joined_lines, to_subset_file, to_subset_file), shell='True')
# for line_to_print in lines_to_sed:
#     subprocess.call("sed -n '{}p' {} >> pc_{}".format(line_to_print, to_subset_file, to_subset_file), shell='True')
