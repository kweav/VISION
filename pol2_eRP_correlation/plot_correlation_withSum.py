#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

'''Usage: date; time ./plot_correlation_withSum.py eRP_pol2.bed eRP_pol2_v.bed'''
gamma = 0.7
num = 11

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

#subselect just the gene_group, eRP value, signal value, num_selected
df_int_sub = df_int.iloc[:,[6,8,15,5]]
#insert the distance correction using ccREstart, ccREend, TSSloc (1-find middle location of ccRE, 2-find abs distance to TSS, 3-scale the distance by ^-gamma)
df_int_sub.insert(4,'dist_correction', pd.DataFrame.abs(((df_int.iloc[:,2]-df_int.iloc[:,1])/2) - df_int.iloc[:,3])**-gamma)
#insert ccRE_name column so that I can find the unique ccREs more easily
df_int_sub.insert(0, 'ccRE_name', df_int[[0, 1, 2]].astype(str).apply(lambda x: "_".join(x), axis=1))


#subselect just the gene_group, eRP value, num_selected
df_v_sub = df_v.iloc[:,[6,8,5]]
#add a signal value of 0
df_v_sub.insert(2,15,0)
#add a distance correction in the same manner as df_int_sub
df_v_sub.insert(4, 'dist_correction', pd.DataFrame.abs(((df_v.iloc[:,2]-df_v.iloc[:,1])/2) - df_v.iloc[:,3])**-gamma)
#insert ccRE_name column so that I can find the unique ccREs more easily
df_v_sub.insert(0, 'ccRE_name', df_v[[0, 1, 2]].astype(str).apply(lambda x: "_".join(x), axis=1))


#combine the two
df_full = pd.concat([df_int_sub, df_v_sub])
