#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

'''Usage: date; time ./plot_correlation.py eRP_pol2.bed eRP_pol2_v.bed'''

def return_cor(x, y):
    cor, pval = stats.pearsonr(x,y)
    coefficients = np.polyfit(x,y,1)
    function_fit = np.poly1d(coefficients)
    x_s = np.linspace(pd.DataFrame.min(x), pd.DataFrame.max(x), 200)
    y_s = function_fit(x_s)
    return(cor, coefficients, x_s, y_s)

def plot_nogg(x,y, fig_name):
    fig, ax = plt.subplots(figsize=(20,5))
    plt.xlim(-40,40)
    ax.set_title('G1E cells correlation between eRP scores and Pol2')
    ax.set_ylabel('Pol2 ChIP Signal Value')
    ax.set_xlabel('eRP score')
    ax.scatter(x, y, alpha=0.5)
    cor, coefficients, x_s, y_s = return_cor(x, y)
    print(cor, coefficients)
    ax.annotate("Pearson R: {0:.3f}".format(cor), (-40, 100))
    ax.plot(x_s, y_s, c='black')
    ax.annotate("{0:.3f}x + {0:.3f}".format(coefficients[0], coefficients[1]), (-40, 75))
    plt.tight_layout()
    fig.savefig(fig_name)
    plt.close(fig)

def plot_gg(x1,y1, x2, y2, x3, y3, x4, y4, fig_name):
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(20,20))
    fig.suptitle('G1E cells correlation between eRP scores and Pol2')
    ax1.set_xlim(-40, 40)
    ax1.set_title('Gene Group 1')
    ax1.scatter(x1, y1, alpha=0.5)
    cor, coefficients, x_s, y_s = return_cor(x1, y1)
    print(cor, coefficients)
    ax1.annotate("Pearson R: {0:.3f}".format(cor), (-40, 100))
    ax1.plot(x_s, y_s, c='black')
    ax1.annotate("{0:.3f}x + {0:.3f}".format(coefficients[0], coefficients[1]), (-40, 75))

    ax2.set_xlim(-40, 40)
    ax2.set_title('Gene Group 2')
    ax2.scatter(x2, y2, alpha=0.5)
    cor, coefficients, x_s, y_s = return_cor(x2, y2)
    print(cor, coefficients)
    ax2.annotate("Pearson R: {0:.3f}".format(cor), (-40, 100))
    ax2.plot(x_s, y_s, c='black')
    ax2.annotate("{0:.3f}x + {0:.3f}".format(coefficients[0], coefficients[1]), (-40, 75))

    ax3.set_xlim(-40, 40)
    ax3.set_title('Gene Group 3')
    ax3.set_ylabel('Pol2 ChIP Signal Value')
    ax3.scatter(x3, y3, alpha=0.5)
    cor, coefficients, x_s, y_s = return_cor(x3, y3)
    print(cor, coefficients)
    ax3.annotate("Pearson R: {0:.3f}".format(cor), (-40, 100))
    ax3.plot(x_s, y_s, c='black')
    ax3.annotate("{0:.3f}x + {0:.3f}".format(coefficients[0], coefficients[1]), (-40, 75))

    ax4.set_xlim(-40,40)
    ax4.set_title('Gene Group 4')
    ax4.scatter(x4, y4, alpha=0.5)
    cor, coefficients, x_s, y_s = return_cor(x4, y4)
    print(cor, coefficients)
    ax4.annotate("Pearson R: {0:.3f}".format(cor), (-40, 100))
    ax4.plot(x_s, y_s, c='black')
    ax4.annotate("{0:.3f}x + {0:.3f}".format(coefficients[0], coefficients[1]), (-40, 75))
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

'''NAIVE correlations'''
'''plot everything that intersects pol2 ChIP regardless of gene group'''
plot_nogg(df_int_sub.iloc[:,1], df_int_sub.iloc[:,2],'naive_correlation_only_int.png')

mask_1 = df_int_sub.iloc[:,0] == 1
mask_2 = df_int_sub.iloc[:,0] == 2
mask_3 = df_int_sub.iloc[:,0] == 3
mask_4 = df_int_sub.iloc[:,0] == 4

'''plot everything that intersects pol2 ChIP subselected by gene group'''
plot_gg(df_int_sub[mask_1].iloc[:,1], df_int_sub[mask_1].iloc[:,2],
        df_int_sub[mask_2].iloc[:,1], df_int_sub[mask_2].iloc[:,2],
        df_int_sub[mask_3].iloc[:,1], df_int_sub[mask_3].iloc[:,2],
        df_int_sub[mask_4].iloc[:,1], df_int_sub[mask_4].iloc[:,2], 'naive_correlation_gene_groups_only_int.png')

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
plot_nogg(df_full.iloc[:,1], df_full.iloc[:,2], 'naive_correlation.png')

'''plot all eRP scores assuming no intersection means signal 0 subselected by gene group '''
plot_gg(df_full[mask_1].iloc[:,1], df_full[mask_1].iloc[:,2],
        df_full[mask_2].iloc[:,1], df_full[mask_2].iloc[:,2],
        df_full[mask_3].iloc[:,1], df_full[mask_3].iloc[:,2],
        df_full[mask_4].iloc[:,1], df_full[mask_4].iloc[:,2], 'naive_correlation_gene_groups.png')

'''Less NAIVE correlations: look at distance''' #what distance function did ABC model use?

'''Less NAIVE correlations: look at number of times selected'''
