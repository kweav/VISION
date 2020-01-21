#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

'''Usage: ./plot_corr_dist_by_chr.py correlations_wrn_{}.txt'''


font = {'size'   : 18}
plt.rc('font', **font)

file_base = sys.argv[1]

chroms = []
for i in range(1,20):
    chroms.append('chr{}'.format(i))
chroms.append('chrX')
colors = ['dodgerblue', 'coral', 'darkslategray', 'purple', 'darkorange', 'olive', 'rosybrown','black', 'khaki','grey', 'tan', 'saddlebrown', 'darkgoldenrod', 'turquoise', 'mediumvioletred', 'cyan', 'navy', 'magenta', 'indigo', 'deeppink']
chrom_to_colors = {}
for color, chrom in zip(colors, chroms):
    chrom_to_colors[chrom] = color

wanted_ticks = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]
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
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7,15))
    ax1.set_title('Distribution of initial Spearmanr\nof expression against ccREs within 1.5Mbp\n{}'.format(chrom))
    ax1.boxplot(chrom_spec_corr)
    ax1.set_yticks(wanted_ticks)
    ax1.set_yticklabels(wanted_ticks)
    ax1.set_ylabel('Spearmanr')
    parts = ax2.violinplot(chrom_spec_corr)
    ax2.set_ylabel('Spearmanr')
    ax2.set_yticks(wanted_ticks)
    ax2.set_yticklabels(wanted_ticks)
    for pc in parts['bodies']:
        pc.set_facecolor(chrom_to_colors[chrom])
        pc.set_edgecolor('black')
        pc.set_alpha(1)
    ax3.hist(chrom_spec_corr, bins=75, color=chrom_to_colors[chrom])
    ax3.set_xlabel('Spearmanr')
    plt.tight_layout()
    fig.savefig('dist_correlation_{}.pdf'.format(chrom), format='pdf')
    fig.savefig('dist_correlation_{}.png'.format(chrom), format='png')
    plt.close(fig)
    #total_corrs += corrs.shape[0]
    #num_nan_total += np.sum(np.isnan(corrs))

correlations = np.array(correlations)
subset = ~np.isnan(correlations)
correlations = correlations[subset]
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7,15))
ax1.set_title('Distribution of initial Spearmanr\nof expression against ccREs within 1.5Mbp\nall chromosomes')
#boxplot
ax1.boxplot(correlations)
ax1.set_yticks(wanted_ticks)
ax1.set_yticklabels(wanted_ticks)
ax1.set_ylabel('Spearmanr')
#violinplot
ax2.violinplot(correlations)
ax2.set_ylabel('Spearmanr')
ax2.set_yticks(wanted_ticks)
ax2.set_yticklabels(wanted_ticks)
#histogram
ax3.hist(correlations, bins=75)
ax3.set_xlabel('Spearmanr')
plt.tight_layout()
fig.savefig('dist_correlation_all.pdf', format='pdf')
fig.savefig('dist_correlation_all.png', format='png')
plt.close(fig)
#print(num_nan_total) #>1159389
#print(total_corrs) #>11965860
