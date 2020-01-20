#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

'''Usage: ./plot_corr_dist_by_chr.py correlations_wrn_{}.txt'''


font = {'size'   : 35}
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

correlations = []
#total_corrs = 0
#num_nan_total = 0
for chrom in chroms:
    chrom_spec_corr = []
    corrs = open(file_base.format(chrom)).readlines()
    for x in corrs:
        correlations.append(float(x))
        chrom_spec_corr.append(float(x))
    chrom_spec_corr = np.array(chrom_spec_corr)
    subset = ~np.isnan(chrom_spec_corr)
    chrom_spec_corr = chrom_spec_corr[subset]
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(30,60))
    fig.suptitle('Distribution of initial Spearmanr correlations\nof expression against ccREs within 1.5Mbp\n{} (subset to include no fully state 0)'.format(chrom))
    ax1.boxplot(chrom_spec_corr)
    ax1.set_ylabel('Spearmanr')
    parts = ax2.violinplot(chrom_spec_corr)
    ax2.set_ylabel('Spearmanr')
    for pc in parts['bodies']:
        pc.set_facecolor(chrom_to_colors[chrom])
        pc.set_edgecolor('black')
        pc.set_alpha(1)
    ax3.hist(chrom_spec_corr, bins=75, color=chrom_to_colors[chrom])
    ax3.set_xlabel('Spearmanr')
    fig.savefig('dist_correlation_{}.png'.format(chrom))
    plt.close(fig)
    #total_corrs += corrs.shape[0]
    #num_nan_total += np.sum(np.isnan(corrs))

correlations = np.array(correlations)
subset = ~np.isnan(correlations)
correlations = correlations[subset]
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(30,60))
fig.suptitle('Distribution of initial Spearmanr correlations\nof expression against ccREs within 1.5Mbp\nall chromosomes (subset to include no fully state 0)')
#boxplot
ax1.boxplot(correlations)
ax1.set_ylabel('Spearmanr')
#violinplot
ax2.violinplot(correlations)
ax2.set_ylabel('Spearmanr')
#histogram
ax3.hist(correlations, bins=75)
ax3.set_xlabel('Spearmanr')
fig.savefig('dist_correlation_all.png')
plt.close(fig)
#print(num_nan_total) #>1159389
#print(total_corrs) #>11965860
