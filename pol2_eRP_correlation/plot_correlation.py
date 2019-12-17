#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

'''Usage: date; time ./plot_correlation.py eRP_pol2.bed eRP_pol2_v.bed'''
gamma = 1
num = 4

def annotate_cor(x, y, ax, an):
    cor, pval = stats.pearsonr(x,y)
    coefficients = np.polyfit(x,y,1)
    print(cor, coefficients)
    function_fit = np.poly1d(coefficients)
    x_s = np.linspace(pd.DataFrame.min(x), pd.DataFrame.max(x), 200)
    y_s = function_fit(x_s)
    ax.annotate("Pearson R: {0:.3f}".format(cor), (-an, 100))
    ax.plot(x_s, y_s, c='black')
    ax.annotate("{0:.3f}x + {0:.3f}".format(coefficients[0], coefficients[1]), (-an, 75))



def plot_nogg(x,y, fig_name, an, xlim = False):
    fig, ax = plt.subplots(figsize=(20,5))
    if xlim:
        plt.xlim(-an,an)
    ax.set_title('G1E cells correlation between eRP scores and Pol2')
    ax.set_ylabel('Pol2 ChIP Signal Value')
    ax.set_xlabel('eRP score')
    ax.scatter(x, y, alpha=0.5)
    annotate_cor(x, y, ax, an)
    plt.tight_layout()
    fig.savefig(fig_name)
    plt.close(fig)

def plot_gg(x1,y1, x2, y2, x3, y3, x4, y4, fig_name, an, xlim = False):
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(20,20))
    fig.suptitle('G1E cells correlation between eRP scores and Pol2')
    if xlim:
        ax1.set_xlim(-an, an)
        ax2.set_xlim(-an, an)
        ax3.set_xlim(-an, an)
        ax4.set_xlim(-an, an)
    ax1.set_title('Gene Group 1')
    ax1.scatter(x1, y1, alpha=0.5)
    annotate_cor(x1, y1, ax1, an)


    ax2.set_title('Gene Group 2')
    ax2.scatter(x2, y2, alpha=0.5)
    annotate_cor(x2, y2, ax2, an)


    ax3.set_title('Gene Group 3')
    ax3.set_ylabel('Pol2 ChIP Signal Value')
    ax3.scatter(x3, y3, alpha=0.5)
    annotate_cor(x3, y3, ax3, an)


    ax4.set_title('Gene Group 4')
    ax4.scatter(x4, y4, alpha=0.5)
    annotate_cor(x4, y4, ax4, an)

    ax4.set_xlabel('eRP score')

    fig.savefig(fig_name)
    plt.tight_layout()
    plt.close(fig)



font = {'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 22}

plt.rc('font', **font)

intersect, non_intersect = sys.argv[1:3]

#gene groups in column 6, eRP_G1E in column 8, signal value in column 15
df_int = pd.read_csv(intersect, sep='\t', header=None)
df_int = df_int.dropna()

#gene groups in column 6, eRP_G1E in column 8, no signal value - assume 0?
df_v = pd.read_csv(non_intersect, sep='\t', header=None)
df_v = df_v.dropna()

#subselect just the gene_group, eRP value, and signal value
df_int_sub = df_int.iloc[:,[6,8,15]]

'''NAIVE CORRELATIONS'''
'''plot everything that intersects pol2 ChIP regardless of gene group'''
#plot_nogg(df_int_sub.iloc[:,1], df_int_sub.iloc[:,2],'naive_correlation_only_int.png', 40, xlim=True)

mask_1 = df_int_sub.iloc[:,0] == 1
mask_2 = df_int_sub.iloc[:,0] == 2
mask_3 = df_int_sub.iloc[:,0] == 3
mask_4 = df_int_sub.iloc[:,0] == 4

'''plot everything that intersects pol2 ChIP subselected by gene group'''
#plot_gg(df_int_sub[mask_1].iloc[:,1], df_int_sub[mask_1].iloc[:,2],
#        df_int_sub[mask_2].iloc[:,1], df_int_sub[mask_2].iloc[:,2],
#        df_int_sub[mask_3].iloc[:,1], df_int_sub[mask_3].iloc[:,2],
#        df_int_sub[mask_4].iloc[:,1], df_int_sub[mask_4].iloc[:,2], 'naive_correlation_gene_groups_only_int.png', 40, xlim=True)

#subselect just the gene_group and eRP value
df_v_sub = df_v.iloc[:,[6,8]]
#add a signal value of 0
df_v_sub.insert(2,15,0)

#combine the two
df_full = pd.concat([df_int_sub, df_v_sub]) #There are NaN values - why??

mask_1 = df_full.iloc[:,0] == 1
mask_2 = df_full.iloc[:,0] == 2
mask_3 = df_full.iloc[:,0] == 3
mask_4 = df_full.iloc[:,0] == 4

'''plot all eRP scores assuming no intersection means signal 0 regardless of gene group '''
#plot_nogg(df_full.iloc[:,1], df_full.iloc[:,2], 'naive_correlation.png', 40, xlim=True)

'''plot all eRP scores assuming no intersection means signal 0 subselected by gene group '''
#plot_gg(df_full[mask_1].iloc[:,1], df_full[mask_1].iloc[:,2],
#        df_full[mask_2].iloc[:,1], df_full[mask_2].iloc[:,2],
#        df_full[mask_3].iloc[:,1], df_full[mask_3].iloc[:,2],
#        df_full[mask_4].iloc[:,1], df_full[mask_4].iloc[:,2], 'naive_correlation_gene_groups.png', 40, xlim=True)

'''LESS NAIVE CORRELATIONS: LOOK AT NUMBER SELECTED'''


'''LESS NAIVE CORRELATIONS: LOOK AT DISTANCE CORRECTION''' #what distance function did ABC model use?
'''ABC model paper talks about how contact could be approximated by distance^-gamma where gamma was 0.7 or 1 for the Extrusion globule model or Fractal globule module respectively. These had similar AUPRCs'''
#so here I'm finding the middle location of the ccRE
#then I'm finding the abs distance from the middle location to the TSS
#finally I'm scaling this distance by ^-gamma like the ABC paper

df_int_sub_dist = df_int_sub.copy()
df_int_sub_dist.insert(3, 'dist_correction', pd.DataFrame.abs(((df_int.iloc[:,2]-df_int.iloc[:,1])/2) - df_int.iloc[:,3])**-gamma)
df_v_sub_dist = df_v_sub.copy()
df_v_sub_dist.insert(3, 'dist_correction', pd.DataFrame.abs(((df_v.iloc[:,2]-df_v.iloc[:,1])/2) - df_v.iloc[:,3])**-gamma)
df_full_dist = pd.concat([df_int_sub_dist, df_v_sub_dist])
'''plot everything that intersects pol2 ChIP regardless of gene group'''
#plot_nogg(df_int_sub_dist.iloc[:,1]*df_int_sub_dist.iloc[:,3], df_int_sub_dist.iloc[:,2], 'abc_dist_correlation_only_int.png', 0.002)

#plot_nogg(df_int_sub_dist.iloc[:,1]/df_int_sub_dist.iloc[:,3], df_int_sub_dist.iloc[:,2], 'abc_dist_div_correlation_only_int.png', 0.002)


mask_1 = df_int_sub_dist.iloc[:,0] == 1
mask_2 = df_int_sub_dist.iloc[:,0] == 2
mask_3 = df_int_sub_dist.iloc[:,0] == 3
mask_4 = df_int_sub_dist.iloc[:,0] == 4

#plot_gg(df_int_sub_dist[mask_1].iloc[:,1]*df_int_sub_dist[mask_1].iloc[:,3], df_int_sub_dist[mask_1].iloc[:,2],
#        df_int_sub_dist[mask_2].iloc[:,1]*df_int_sub_dist[mask_2].iloc[:,3], df_int_sub_dist[mask_2].iloc[:,2],
#        df_int_sub_dist[mask_3].iloc[:,1]*df_int_sub_dist[mask_3].iloc[:,3], df_int_sub_dist[mask_3].iloc[:,2],
#        df_int_sub_dist[mask_4].iloc[:,1]*df_int_sub_dist[mask_4].iloc[:,3], df_int_sub_dist[mask_4].iloc[:,2], 'abc_dist_correlation_gene_groups_only_int.png', 0.002)

#plot_gg(df_int_sub_dist[mask_1].iloc[:,1]/df_int_sub_dist[mask_1].iloc[:,3], df_int_sub_dist[mask_1].iloc[:,2],
#        df_int_sub_dist[mask_2].iloc[:,1]/df_int_sub_dist[mask_2].iloc[:,3], df_int_sub_dist[mask_2].iloc[:,2],
#        df_int_sub_dist[mask_3].iloc[:,1]/df_int_sub_dist[mask_3].iloc[:,3], df_int_sub_dist[mask_3].iloc[:,2],
#        df_int_sub_dist[mask_4].iloc[:,1]/df_int_sub_dist[mask_4].iloc[:,3], df_int_sub_dist[mask_4].iloc[:,2], 'abc_dist_div_correlation_gene_groups_only_int.png', 0.002)

'''plot all eRP scores assuming no intersection means signal 0 regardless of gene group '''
#plot_nogg(df_full_dist.iloc[:,1]*df_full_dist.iloc[:,3], df_full_dist.iloc[:,2], 'abc_dist_correlation.png', 0.002)

#plot_nogg(df_full_dist.iloc[:,1]/df_full_dist.iloc[:,3], df_full_dist.iloc[:,2], 'abc_dist_div_correlation.png', 0.002)

mask_1 = df_full_dist.iloc[:,0] == 1
mask_2 = df_full_dist.iloc[:,0] == 2
mask_3 = df_full_dist.iloc[:,0] == 3
mask_4 = df_full_dist.iloc[:,0] == 4

#plot_gg(df_full_dist[mask_1].iloc[:,1]*df_full_dist[mask_1].iloc[:,3], df_full_dist[mask_1].iloc[:,2],
#        df_full_dist[mask_2].iloc[:,1]*df_full_dist[mask_2].iloc[:,3], df_full_dist[mask_2].iloc[:,2],
#        df_full_dist[mask_3].iloc[:,1]*df_full_dist[mask_3].iloc[:,3], df_full_dist[mask_3].iloc[:,2],
#        df_full_dist[mask_4].iloc[:,1]*df_full_dist[mask_4].iloc[:,3], df_full_dist[mask_4].iloc[:,2], 'abc_dist_correlation_gene_groups.png', 0.002)

#plot_gg(df_full_dist[mask_1].iloc[:,1]/df_full_dist[mask_1].iloc[:,3], df_full_dist[mask_1].iloc[:,2],
#        df_full_dist[mask_2].iloc[:,1]/df_full_dist[mask_2].iloc[:,3], df_full_dist[mask_2].iloc[:,2],
#        df_full_dist[mask_3].iloc[:,1]/df_full_dist[mask_3].iloc[:,3], df_full_dist[mask_3].iloc[:,2],
#        df_full_dist[mask_4].iloc[:,1]/df_full_dist[mask_4].iloc[:,3], df_full_dist[mask_4].iloc[:,2], 'abc_dist_div_correlation_gene_groups.png', 0.002)


#Mike's concern: What about a cCRE connected to multiple genes at the same time. That's a good question. How do I handle that? What should I be averaging?

'''LESS NAIVE CORRELATIONS: LOOK AT NUMBER OF TIMES SELECTED'''
