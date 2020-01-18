#!/usr/bin/env python3

'''Usage: ./num_in_window.py file_base_{}.txt png'''

import sys
import numpy as np
import matplotlib.pyplot as plt

font = {'size'   : 22}

plt.rc('font', **font)

file_base = sys.argv[1]

chroms = []
for i in range(1,19):
    chroms.append('chr{}'.format(i))
chroms.append('chrX')
colors = ['dodgerblue', 'coral', 'darkslategray', 'purple', 'darkorange', 'olive', 'rosybrown','black', 'khaki','grey', 'saddlebrown', 'tan', 'darkgoldenrod', 'turquoise', 'mediumvioletred', 'cyan', 'navy', 'indigo', 'magenta', 'deeppink']

chrom_to_colors = {}
for color, chrom in zip(colors, chroms):
    chrom_to_colors[chrom] = color

num_in_window = []
colors_to_plot = []
for chrom in chroms:
    for i, line in enumerate(open(file_base.format(chrom))):
        num = int(line.strip('\r\n').split('\t')[1])
        if num == 0:
            print(chrom, i, flush=True)
        num_in_window.append(num)
        colors_to_plot.append(chrom_to_colors[chrom])

numTSS = np.arange(len(num_in_window))

fig, ax = plt.subplots(figsize=(20,20))
fig.suptitle('#ccREs within 1.5Mbp of each TSS (colored by chr)')
ax.scatter(numTSS, num_in_window, c=colors_to_plot)
ax.set_xlabel('TSSn')
ax.set_ylabel('#ccREs')
fig.savefig(sys.argv[2])
plt.close(fig)
