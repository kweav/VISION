#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

font = {'size'   : 20}
plt.rc('font', **font)

npzfile_base = sys.argv[1]

chroms = []
for value in range(1,20):
    chroms.append('chr{}'.format(value))
chroms.append('chrX')

num_pairings = np.zeros((len(chroms), 13), dtype=np.object)
for i, chrom in enumerate(chroms):
    npzfile = np.load(npzfile_base.format(chrom))#, allow_pickle = True)
    num_selected = npzfile['num_selected']
    num_pairings[i, 0] = np.sum(~np.isnan(num_selected), axis=1)
    for j in range(1,13):
        num_pairings[i, j] = np.sum(num_selected == j, axis=1)

datasets_to_plot = []
labels = []
for k in range(13):
    if k == 0:
        labels.append('Filter\nSet')
    else:
        labels.append('{}'.format(k))
    datasets_to_plot.append(np.hstack(num_pairings[:,k]))

fig, ax = plt.subplots(figsize=(30,7))
ax.set_title('Predicted Enhancer-Promoter Interaction')
parts = ax.violinplot(datasets_to_plot, positions=np.arange(len(datasets_to_plot)))
ax.boxplot(datasets_to_plot, positions=np.arange(len(datasets_to_plot)))
ax.set_ylabel("# of enhancers paired with each promoter")
ax.set_xlabel("# of training iterations pairing selected")
ax.set_xticks(np.arange(len(datasets_to_plot)))
ax.set_xticklabels(labels)
for pc in parts['bodies']:
    pc.set_facecolor('C0')
    pc.set_edgecolor('C0')
    pc.set_alpha(1)
plt.tight_layout()
fig.savefig("allChr_numpaired_per_gene.png")
plt.close(fig)

fig, ax = plt.subplots(figsize=(30,7))
ax.set_title('Predicted Enhancer-Promoter Interaction')
parts = ax.violinplot(datasets_to_plot[1:], positions=np.arange(len(datasets_to_plot[1:])))
ax.boxplot(datasets_to_plot[1:], positions=np.arange(len(datasets_to_plot[1:])))
ax.set_ylabel("# of enhancers paired with each promoter")
ax.set_xlabel("# of training iterations pairing selected")
ax.set_xticks(np.arange(len(datasets_to_plot[1:])))
ax.set_xticklabels(labels[1:])
for pc in parts['bodies']:
    pc.set_facecolor('C0')
    pc.set_edgecolor('C0')
    pc.set_alpha(1)
plt.tight_layout()
fig.savefig("allChr_numpaired_per_gene_nfs.png")
plt.close(fig)
