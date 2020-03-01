#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle

font = {'size'   : 20}
plt.rc('font', **font)

datasets = {"full":
                    {'r1: 0-11': [],
                     'r2: 0-11': []},
            "goe": {'r1: 0-11': [],
                    'r2: 0-11':[]}}
ordered_keys = ['r1: 0-11', 'r2: 0-11']
ordered_goe_keys = ['r1: 0-11', 'r2: 0-11']

for i in range(13):
    datasets['full']['r1: {}'.format(i)] = []
    datasets['full']['r2: {}'.format(i)] = []
    ordered_keys.append('r1: {}'.format(i))
    ordered_keys.append('r2: {}'.format(i))
    datasets['goe']['r1: >={}'.format(i)] = []
    datasets['goe']['r2: >={}'.format(i)] = []
    ordered_goe_keys.append('r1: >={}'.format(i))
    ordered_goe_keys.append('r2: >={}'.format(i))

file_base = 'subset_{}.pickle'
for i in range(1, 28):
    with open(file_base.format(i), 'rb') as f:
        pickled = pickle.load(f, encoding='latin1')
        for label, dataset_list in zip(pickled['labels'], pickled['datasets']):
            for value in dataset_list:
                datasets['full'][label].append(float(value))
        for label, dataset_list in zip(pickled["labels_goe"], pickled["datasets_goe"]):
            for value in dataset_list:
                datasets['goe'][label].append(float(value))
with open(file_base.format('final'), 'rb') as f:
    pickled = pickle.load(f, encoding='latin1')
    for label, dataset_list in zip(pickled['labels'], pickled['datasets']):
        for value in dataset_list:
            datasets['full'][label].append(float(value))
    for label, dataset_list in zip(pickled["labels_goe"], pickled["datasets_goe"]):
        for value in dataset_list:
            datasets['goe'][label].append(float(value))


datasets_full = []
labels_full = []
for key in ordered_keys:
    if 'r2' not in key:
        full_to_subset = np.array(datasets['full'][key])
        full = full_to_subset[~np.isnan(full_to_subset)]
        datasets_full.append(full)

        if key == 'r1: 0-11':
            labels_full.append('Filter\nSet')
        else:
            labels_full.append(key.replace('r1: ', ''))


labels_full = np.array(labels_full)
datasets_full = np.array(datasets_full)

# fig, ax = plt.subplots(figsize=(30,20))
# ax.set_title('All Chr (random subset) G1e Interaction Frequencies')
# ax.boxplot(datasets_full, showfliers = False)
# ax.set_ylabel("Enrichment Interaction Frequency (observed/expected)")
# ax.set_xlabel("Number of times selected for inclusion")
# #ax.set_xticks(np.arange(len(datasets_full)))
# ax.set_xticklabels(labels_full)
# fig.savefig("box_subset_G1e_int_freq_num_nofliers.png")
# plt.close(fig)
#
# fig, ax = plt.subplots(figsize=(30,20))
# ax.set_title('Enhancer-Promoter Interaction')
#
# ax.set_ylabel("Interaction Frequency Enrichment")
# ax.set_xlabel("number of times selected")
# #ax.set_xticks(np.arange(len(datasets_full)))
# ax.set_xticklabels(labels_full)
# fig.savefig("box_subset_G1e_int_freq_num.png")
# plt.close(fig)
#
# fig, ax = plt.subplots(figsize=(30,20))
# ax.set_title('All Chr (random subset) G1e Interaction Frequencies')
# ax.boxplot(datasets_goe)
# ax.set_ylabel("Enrichment Interaction Frequency (observed/expected)")
# ax.set_xlabel("Number of times selected for inclusion")
# #ax.set_xticks(np.arange(len(datasets_goe)))
# ax.set_xticklabels(labels_goe)
# fig.savefig("box_subset_G1e_int_freq_goe.png")
# plt.close(fig)
#
# fig, ax = plt.subplots(figsize=(30,20))
# ax.set_title('All Chr (random subset) G1e Interaction Frequencies')
# ax.boxplot(datasets_goe, showfliers=False)
# ax.set_ylabel("Enrichment Interaction Frequency (observed/expected)")
# ax.set_xlabel("Number of times selected for inclusion")
# #ax.set_xticks(np.arange(len(datasets_goe)))
# ax.set_xticklabels(labels_goe)
# fig.savefig("box_subset_G1e_int_freq_goe_nofliers.png")
# plt.close(fig)
#
# fig, ax = plt.subplots(figsize=(30,20))
# ax.set_title('All Chr (random subset) G1e Interaction Frequencies')
# ax.violinplot(datasets_full, showextrema = False)
# ax.set_ylabel("Enrichment Interaction Frequency (observed/expected)")
# ax.set_xlabel("Number of times selected for inclusion")
# #ax.set_xticks(np.arange(len(datasets_full)))
# ax.set_xticklabels(labels_full)
# fig.savefig("vi_subset_G1e_int_freq_num_nofliers.png")
# plt.close(fig)

fig, ax = plt.subplots(figsize=(30,7))
ax.set_title('Predicted Enhancer-Promoter Interaction')
parts = ax.violinplot(datasets_full, positions=np.arange(len(datasets_full)))
ax.boxplot(datasets_full, positions=np.arange(len(datasets_full)))
ax.set_ylabel("Interaction Frequency Enrichment")
ax.set_xlabel("number of training iterations pairing selected")
ax.set_xticks(np.arange(len(datasets_full)))
ax.set_xticklabels(labels_full)
for pc in parts['bodies']:
    pc.set_facecolor('C0')
    pc.set_edgecolor('C0')
    pc.set_alpha(1)
plt.tight_layout()
fig.savefig("both_subset_G1e_int_freq_num.png")
plt.close(fig)

# fig, ax = plt.subplots(figsize=(30,20))
# ax.set_title('All Chr (random subset) G1e Interaction Frequencies')
# ax.violinplot(datasets_goe)
# ax.set_ylabel("Enrichment Interaction Frequency (observed/expected)")
# ax.set_xlabel("Number of times selected for inclusion")
# #ax.set_xticks(np.arange(len(datasets_goe)))
# ax.set_xticklabels(labels_goe)
# fig.savefig("vi_subset_G1e_int_freq_goe.png")
# plt.close(fig)
#
# fig, ax = plt.subplots(figsize=(30,20))
# ax.set_title('All Chr (random subset) G1e Interaction Frequencies')
# ax.violinplot(datasets_goe, showextrema=False)
# ax.set_ylabel("Enrichment Interaction Frequency (observed/expected)")
# ax.set_xlabel("Number of times selected for inclusion")
# #ax.set_xticks(np.arange(len(datasets_goe)))
# ax.set_xticklabels(labels_goe)
# fig.savefig("vi_subset_G1e_int_freq_goe_nofliers.png")
# plt.close(fig)
