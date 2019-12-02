#!/usr/bin/env python3

'''Usage: ./plot_gene_group_beta.py tss_eRP_2kb.with0.txt dist_eRP_2kb.with0.txt'''

import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

tss,dist = sys.argv[1:]
df_tss = pd.read_csv(tss, sep="\t")
tss_header = list(df_tss.columns.values)
df_tss.columns = ["TSS_" + x for x in tss_header]
df_dist = pd.read_csv(dist, sep='\t')
dist_header = list(df_dist.columns.values)
df_dist.columns = ['cCRE_'+ x for x in dist_header]

df_concat = pd.concat([df_tss, df_dist], axis=1, sort=False)
#print(df_concat)

reorder = [7,4,3,25,2,16,0,1,14,20,13,11,24,23,6,15,12,21,26,10,19,8,18,5,22,9,17]
df_concat_reindex = df_concat.reindex(reorder)
print(df_concat_reindex)

for i in range(8):
    scalar = df_concat_reindex.iloc[6,i]
    df_concat_reindex.iloc[:,i] -= scalar

print(df_concat_reindex)

ax = sns.heatmap(df_concat_reindex, vmin=-40, vmax=40, linewidths=0.25, linecolor='gray', cmap='bwr', yticklabels=1, square=True)#, center=df_concat_reindex.iloc[6,0])
ax.set_yticklabels(list(df_concat_reindex.index), rotation=0)
ax.yaxis.tick_right()
ax.set_xticklabels([x.replace('_gene_c','') for x in df_concat_reindex.columns.values], rotation=270)
plt.tight_layout()
fig = ax.get_figure()
fig.savefig('beta_heatmap_sub0.png')
