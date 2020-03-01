#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

font = {'size'   : 20}
plt.rc('font', **font)

npzfile_base = sys.argv[1]

cellTypes = ['B', 'CD4', 'CD8', 'G1E', 'HPC7', 'LSK', 'Neu']

chroms = []
for value in range(1,20):
    chroms.append('chr{}'.format(value))
chroms.append('chrX')

avg_datasets_to_plot = []
labels_avg = []
num_pairings_avg = np.zeros((len(chroms), 13), dtype=np.object)
for i, chrom in enumerate(chroms):
    npzfile = np.load(npzfile_base.format(chrom))#, allow_pickle = True)
    num_selected = npzfile['num_selected']
    int_freq = npzfile['int_freq']
    num_pairings_avg[i, 0] = np.mean(int_freq[~np.isnan(num_selected)], axis=1)
    for j in range(1,13):
        num_pairings_avg[i,j] = np.mean(int_freq[num_selected==j], axis=1)
for k in range(13):
    if k  == 0:
        labels_avg.append('Filter\nSet')
    else:
        labels_avg.append('{}'.format(k))
    avg_datasets_to_plot.append(np.hstack(num_pairings_avg[:,k]))

fig, ax = plt.subplots(figsize=(30,7))
ax.set_title('Predicted Enhancer-Promoter Interaction')
parts = ax.violinplot(avg_datasets_to_plot[1:], positions=np.arange(len(avg_datasets_to_plot[1:])))
ax.boxplot(avg_datasets_to_plot[1:], positions=np.arange(len(avg_datasets_to_plot[1:])))
ax.set_ylabel("Interaction Frequency Enrichment")
ax.set_xlabel("# of training iterations pairing selected")
ax.set_xticks(np.arange(len(avg_datasets_to_plot[1:])))
ax.set_xticklabels(labels_avg[1:])
for pc in parts['bodies']:
    pc.set_facecolor('C0')
    pc.set_edgecolor('C0')
    pc.set_alpha(1)
plt.tight_layout()
fig.savefig("allChr_int_freq_nfs_avgCT.png")
plt.close(fig)

fig, ax = plt.subplots(figsize=(30,7))
ax.set_title('Predicted Enhancer-Promoter Interaction')
parts = ax.violinplot(avg_datasets_to_plot, positions=np.arange(len(avg_datasets_to_plot)))
ax.boxplot(avg_datasets_to_plot, positions=np.arange(len(avg_datasets_to_plot)))
ax.set_ylabel("Interaction Frequency Enrichment")
ax.set_xlabel("# of training iterations pairing selected")
ax.set_xticks(np.arange(len(avg_datasets_to_plot)))
ax.set_xticklabels(labels_avg)
for pc in parts['bodies']:
    pc.set_facecolor('C0')
    pc.set_edgecolor('C0')
    pc.set_alpha(1)
plt.tight_layout()
fig.savefig("allChr_int_freq_avgCT.png")
plt.close(fig)

# for n, CT in enumerate(cellTypes):
#     num_pairings = np.zeros((len(chroms), 13), dtype=np.object)
#     for i, chrom in enumerate(chroms):
#         npzfile = np.load(npzfile_base.format(chrom))#, allow_pickle = True)
#         num_selected = npzfile['num_selected']
#         int_freq = npzfile['int_freq']
#         num_pairings[i, 0] = int_freq[:,:,n][~np.isnan(num_selected)]
#         for j in range(1,13):
#             num_pairings[i, j] = int_freq[:,:,n][num_selected == j]
#     datasets_to_plot = []
#     labels = []
#     for k in range(13):
#         if k == 0:
#             labels.append('Filter\nSet')
#         else:
#             labels.append('{}'.format(k))
#         datasets_to_plot.append(np.hstack(num_pairings[:,k]))
#
#
#     fig, ax = plt.subplots(figsize=(30,7))
#     ax.set_title('Predicted Enhancer-Promoter Interaction - {} Cells'.format(CT))
#     parts = ax.violinplot(datasets_to_plot[1:], positions=np.arange(len(datasets_to_plot[1:])))
#     ax.boxplot(datasets_to_plot[1:], positions=np.arange(len(datasets_to_plot[1:])))
#     ax.set_ylabel("Interaction Frequency Enrichment")
#     ax.set_xlabel("# of training iterations pairing selected")
#     ax.set_xticks(np.arange(len(datasets_to_plot[1:])))
#     ax.set_xticklabels(labels[1:])
#     for pc in parts['bodies']:
#         pc.set_facecolor('C0')
#         pc.set_edgecolor('C0')
#         pc.set_alpha(1)
#     plt.tight_layout()
#     fig.savefig("allChr_int_freq_nfs_{}.png".format(CT))
#     plt.close(fig)
#
#     fig, ax = plt.subplots(figsize=(30,7))
#     ax.set_title('Predicted Enhancer-Promoter Interaction - {} Cells'.format(CT))
#     parts = ax.violinplot(datasets_to_plot, positions=np.arange(len(datasets_to_plot)))
#     ax.boxplot(datasets_to_plot, positions=np.arange(len(datasets_to_plot)))
#     ax.set_ylabel("Interaction Frequency Enrichment")
#     ax.set_xlabel("# of training iterations pairing selected")
#     ax.set_xticks(np.arange(len(datasets_to_plot)))
#     ax.set_xticklabels(labels)
#     for pc in parts['bodies']:
#         pc.set_facecolor('C0')
#         pc.set_edgecolor('C0')
#         pc.set_alpha(1)
#     plt.tight_layout()
#     fig.savefig("allChr_int_freq_{}.png".format(CT))
#     plt.close(fig)
