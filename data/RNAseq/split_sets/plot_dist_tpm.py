#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import numpy as np

values_train = []
for i, line in enumerate(open(sys.argv[1])):
    if i !=0:
        fields=line.strip('\r\n').split()[6:]
        for value in fields:
            values_train.append(float(value))

train_max = max(values_train)

values_test = []
for i, line in enumerate(open(sys.argv[2])):
    if i !=0:
        fields=line.strip('\r\n').split()[6:]
        for value in fields:
            values_test.append(float(value))

test_max = max(values_test)

values_ref = []
for i, line in enumerate(open(sys.argv[3])):
    if i !=0:
        fields=line.strip('\r\n').split()[6:]
        for value in fields:
            values_ref.append(float(value))

ref_max = max(values_ref)

max_val = np.ceil(max(train_max, test_max, ref_max))

def plot_vals(to_plot, title, figname, max_val=max_val):
    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.hist(to_plot, bins=100)
    ax.set_xlim(-10, max_val)
    ax.set_xlabel('log2(TPM+0.001)')
    ax.set_ylabel('number of genes')
    plt.tight_layout()
    fig.savefig(figname)
    plt.close(fig)


plot_vals(values_train, 'Expression of Training Set of Protein Coding Genes', 'train_dist_tpm_pc.png')
plot_vals(values_test, 'Expression of Testing Set of Protein Coding Genes', 'test_dist_tpm_pc.png')
plot_vals(values_ref, 'Expression of Reference Set of Protein Coding Genes', 'ref_dist_tpm_pc.png')
