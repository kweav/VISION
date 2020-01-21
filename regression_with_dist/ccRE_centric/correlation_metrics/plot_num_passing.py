#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import numpy as np

'''Usage: ./plot_num_passing.py correlationShapes_{}_{}.txt noPairings_possible_{}_{}.txt distances_afterPassing_{}_{}.txt bins'''

file_base_num_passing, file_base_not_paired, file_base_distance, bins_arg = sys.argv[1:]

font = {'size'   : 18}
plt.rc('font', **font)

'''Basic Information'''
chroms = []
for i in range(1, 20):
    chroms.append('chr{}'.format(i))
chroms.append('chrX')

colors = ['dodgerblue', 'coral', 'darkslategray', 'purple', 'darkorange', 'olive', 'rosybrown','black', 'khaki','grey', 'tan', 'saddlebrown', 'darkgoldenrod', 'turquoise', 'mediumvioletred', 'cyan', 'navy', 'magenta', 'indigo', 'deeppink']

corrs = [0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.75]

corrs_to_colors = {}
for color, corr in zip(colors, corrs):
    corrs_to_colors[corr] = color

chroms_to_colors  = {}
for color, chrom in zip(colors, chroms):
    chroms_to_colors[chrom] = color

chrom_len = {'chr1' : 195471971,
             'chr2' : 182113224,
             'chrX' : 171031299,
             'chr3' : 160039680,
             'chr4' : 156508116,
             'chr5' : 151834684,
             'chr6' : 149736546,
             'chr7' : 145441459,
             'chr10' : 130694993,
             'chr8' : 129401213,
             'chr14' : 124902244,
             'chr9' : 124595110,
             'chr11' : 122082543,
             'chr13' : 120421639,
             'chr12' : 120129022,
             'chr15' : 104043685,
             'chr16' : 98207768,
             'chr17' : 94987271,
             'chrY' : 91744698,
             'chr18' : 90702639,
             'chr19' : 61431566}

chrom_genes = {'chr1' : 3415,
               'chr2' : 3748,
               'chr3' : 2192,
               'chr4' : 2878,
               'chr5' : 2059,
               'chr6' : 2134,
               'chr7' : 3461,
               'chr8' : 1768,
               'chr9' : 2140,
               'chr10' : 1613,
               'chr11' : 2974,
               'chr12' : 1483,
               'chr13' : 1491,
               'chr14' : 1739,
               'chr15' : 1249,
               'chr16' : 1155,
               'chr17' : 1701,
               'chr18' : 886,
               'chr19' : 1030,
               'chrX' : 2624}


'''Set up data storage'''
info_overload = {}
chrom_to_colors = {}
for corr in corrs:
    info_overload[corr] = {}
    for color, chrom in zip(colors, chroms):
        info_overload[corr][chrom] = {'locs_not_paired': [],
                                      'num_passing': [],
                                      'dists_passing': [],
                                      'num_not_paired': 0}
        chrom_to_colors[chrom] = color
'''Read in data'''
for corr in corrs:
    for chrom in chroms:
        for line in open(file_base_num_passing.format(corr, chrom)):
            info_overload[corr][chrom]['num_passing'].append(int(line.strip('\r\n').split('\t')[2]))
        for line in open(file_base_not_paired.format(corr, chrom)):
            TSS = int(line.strip('\r\n').split('\t')[1].split(':')[1])
            info_overload[corr][chrom]['locs_not_paired'].append(TSS)
            info_overload[corr][chrom]['num_not_paired'] += 1
        info_overload[corr][chrom]['num_not_paired'] = np.array(info_overload[corr][chrom]['num_not_paired'])
        info_overload[corr][chrom]['num_not_paired'] = np.divide(info_overload[corr][chrom]['num_not_paired'], chrom_genes[chrom]) #proportion of genes
        info_overload[corr][chrom]['locs_not_paired'] = np.array(info_overload[corr][chrom]['locs_not_paired'])
        info_overload[corr][chrom]['locs_not_paired'] = np.divide(info_overload[corr][chrom]['locs_not_paired'], chrom_len[chrom]) #scale by length of chromosome to find standardized loc
        for line in open(file_base_distance.format(corr, chrom)):
            info_overload[corr][chrom]['dists_passing'].append(int(float(line.strip('\r\n')))//1000)


'''plot correlation vs num_not_paired and correlation vs distribtion_of_num_passing'''
for chrom in chroms:
    num_not_paired = []
    datasets = []
    for corr in corrs:
        num_not_paired.append(info_overload[corr][chrom]['num_not_paired'])
        datasets.append(info_overload[corr][chrom]['num_passing'])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7,15))
    ax1.set_title('Balancing # ccREs per pairing\n# genes not paired at all\n{}'.format(chrom))
    ax1.plot(np.arange(len(corrs)), num_not_paired, marker='o')
    ax1.set_xticks(np.arange(len(corrs)))
    ax1.set_xticklabels(labels = corrs)
    ax1.set_xlabel('Abs(Correlation) Thresholds')
    ax1.set_ylabel('# genes without initial pairing/total genes on chr')
    ax2.violinplot(datasets, positions=np.arange(len(corrs)))
    ax2.set_xlabel('Abs(Correlation) Thresholds')
    ax2.set_ylabel('# of paired ccREs per TSS')
    ax2.set_xticks(np.arange(len(corrs)))
    ax2.set_xticklabels(labels = corrs)
    fig.savefig('balance_numPerPairing_noPairing_{}.pdf'.format(chrom), format='pdf')
    fig.savefig('balance_numPerPairing_noPairing_{}.png'.format(chrom), format='png')
    plt.tight_layout()
    plt.close(fig)

'''plot distances from TSS that are retained'''
for chrom in chroms:
    fig1, axes1 = plt.subplots(len(corrs), 1)#, figsize=(30,60))
    fig2, ax2 = plt.subplots(figsize=(15,15))
    ax2.set_xlim(-1500, 1500)
    for i, corr in enumerate(corrs):
        if len(info_overload[corr][chrom]['dists_passing']) > 0:
            y, bins, patches = axes1[i].hist(info_overload[corr][chrom]['dists_passing'], bins=int(bins_arg))
            bincenters = 0.5*(bins[1:]+bins[:-1])
            ax2.plot(bincenters, y, color=corrs_to_colors[corr])
        else:
            ax2.plot(np.linspace(-1500,1500, 4000), np.tile(0, 4000))

    ax2.legend(labels=corrs, loc='upper right')
    ax2.set_title('Distance from center of ccRE to initially paired TSS')
    ax2.set_xlabel('Distance from TSS (kbps)')
    ax2.set_ylabel('Number of paired ccREs')
    plt.tight_layout()
    fig2.savefig('distances_represented_{}_bins_{}.pdf'.format(chrom, bins_arg), format='pdf')
    fig2.savefig('distances_represented_{}_bins_{}.png'.format(chrom, bins_arg), format='png')
    plt.close(fig1)
    plt.close(fig2)

'''plot where on the chromosomes the unpaired TSSs are'''
fig, axes = plt.subplots(len(corrs), 1, figsize=(7, 15))
for i, corr in enumerate(corrs):
    for j, chrom in enumerate(chroms, 1):
        axes[i].scatter(info_overload[corr][chrom]['locs_not_paired'], np.tile(j, info_overload[corr][chrom]['locs_not_paired'].shape[0]), c=chroms_to_colors[chrom])
    if i ==0:
        axes[i].set_title('Location of unpaired TSSs\nSpearmanr Threshold: {}'.format(corr))
    else:
        axes[i].set_title('Spearmanr Threshold: {}'.format(corr))
    if i == (len(corrs)-1)//2:
        axes[i].set_ylabel('Chromosome')
    #axes[i].set_yticks(np.arange(1, 21))
    #axes[i].set_yticklabels(chroms)
    if i == (len(corrs) - 1):
        axes[i].set_xlabel('Relative location on chromosome')
    axes[i].set_xlim(0, 1)
    axes[i].set_xticks([0, 0.25, 0.5, 0.75, 1])
    axes[i].set_xticklabels([0, 0.25, 0.5, 0.75, 1])
plt.tight_layout()
fig.savefig('loc_unpaired.pdf', format='pdf')
fig.savefig('loc_unpaired.png', format='png')
plt.close(fig)
