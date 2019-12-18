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
    x_s = np.linspace(np.min(x), np.max(x), 200)
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
df_int_sub.insert(0, 'ccRE_name', df_int[0] + "_" + df_int[1].astype(str) + "_" + df_int[2].astype(str))

unique_ccREs_int = np.array(df_int_sub.ccRE_name.unique())

'''only those that intersect pol2 ChIP'''
summed_eRPs = np.zeros((unique_ccREs_int.shape[0], 3))
pol2_vals = np.zeros(unique_ccREs_int.shape[0])
values_gg = np.full((unique_ccREs_int.shape[0], 4, 4), np.nan) # first dimension is ccRE; second dimension is sum or pol2 val ;third dimension is gg
for i, ccRE in enumerate(unique_ccREs_int):
    '''regardless of gene group'''
    mask = df_int_sub.iloc[:,0] == ccRE
    '''sum of raw'''
    sum_e = pd.DataFrame.sum(df_int_sub[mask].iloc[:,2])
    '''sum of distance corrected'''
    sum_e_dc_p = pd.DataFrame.sum(df_int_sub[mask].iloc[:,2]*df_int_sub[mask].iloc[:,5])
    sum_e_dc_d = pd.DataFrame.sum(df_int_sub[mask].iloc[:,2]/df_int_sub[mask].iloc[:,5])
    summed_eRPs[i,[0,1,2]] = [sum_e, sum_e_dc_p, sum_e_dc_d]
    pol2_val = df_int_sub[mask].iloc[:,3].unique()[0]
    pol2_vals[i] = pol2_val

    # '''taking into account gene group'''
    # mask1 = (df_int_sub.iloc[:,0] == ccRE) & (df_int_sub.iloc[:,1] == 1)
    # mask2 = (df_int_sub.iloc[:,0] == ccRE) & (df_int_sub.iloc[:,1] == 2)
    # mask3 = (df_int_sub.iloc[:,0] == ccRE) & (df_int_sub.iloc[:,1] == 3)
    # mask4 = (df_int_sub.iloc[:,0] == ccRE) & (df_int_sub.iloc[:,1] == 4)
    #
    # sum_e1 = pd.DataFrame.sum(df_int_sub[mask1].iloc[:,2])
    # sum_e_dc_p1 = pd.DataFrame.sum(df_int_sub[mask1].iloc[:,2]*df_int_sub[mask1].iloc[:,5])
    # sum_e_dc_d1 = pd.DataFrame.sum(df_int_sub[mask1].iloc[:,2]/df_int_sub[mask1].iloc[:,5])
    # pol2_val = df_int_sub[mask1].iloc[:,3].unique()[0]
    # values_gg[i, [0,1,2,3], 0] = [sum_e1, sum_e_dc_p1, sum_e_dc_d1, pol2_val]
    #
    # sum_e2 = pd.DataFrame.sum(df_int_sub[mask2].iloc[:,2])
    # sum_e_dc_p2 = pd.DataFrame.sum(df_int_sub[mask2].iloc[:,2]*df_int_sub[mask2].iloc[:,5])
    # sum_e_dc_d2 = pd.DataFrame.sum(df_int_sub[mask2].iloc[:,2]/df_int_sub[mask2].iloc[:,5])
    # pol2_val = df_int_sub[mask2].iloc[:,3].unique()[0]
    # values_gg[i, [0,1,2,3], 1] = [sum_e2, sum_e_dc_p2, sum_e_dc_d2, pol2_val]
    #
    # sum_e3 = pd.DataFrame.sum(df_int_sub[mask3].iloc[:,2])
    # sum_e_dc_p3 = pd.DataFrame.sum(df_int_sub[mask3].iloc[:,2]*df_int_sub[mask3].iloc[:,5])
    # sum_e_dc_d3 = pd.DataFrame.sum(df_int_sub[mask3].iloc[:,2]/df_int_sub[mask3].iloc[:,5])
    # pol2_val = df_int_sub[mask3].iloc[:,3].unique()[0]
    # values_gg[i, [0,1,2,3], 2] = [sum_e3, sum_e_dc_p3, sum_e_dc_d3, pol2_val]
    #
    # sum_e4 = pd.DataFrame.sum(df_int_sub[mask4].iloc[:,2])
    # sum_e_dc_p4 = pd.DataFrame.sum(df_int_sub[mask4].iloc[:,2]*df_int_sub[mask4].iloc[:,5])
    # sum_e_dc_d4 = pd.DataFrame.sum(df_int_sub[mask4].iloc[:,2]/df_int_sub[mask4].iloc[:,5])
    # pol2_val = df_int_sub[mask4].iloc[:,3].unique()[0]
    # values_gg[i, [0,1,2,3], 3] = [sum_e4, sum_e_dc_p4, sum_e_dc_d4, pol2_val]

#np.save('only_int_values_gg_gamma{}.npy'.format(gamma), values_gg)

'''Plot summed raw eRP (1st of second dimension of summed_eRPs)
--> only those that intersect pol2 ChIP
--> --> regardless of gene group'''
plot_nogg(summed_eRPs[:,0], pol2_vals, 'sum_raw_onlyint.png', np.max(summed_eRPs[:,0])+5, xlim=True)

'''Plot summed distance corrected (*) eRP (2nd of second dimension of summed_eRPs)
--> only those that intersect pol2 ChIP
--> --> regardless of gene group'''
plot_nogg(summed_eRPs[:,1], pol2_vals, 'sum_dc_p_gamma{}_onlyint.png'.format(gamma), np.max(summed_eRPs[:,1])+5, xlim=True)

'''Plot summed distance corrected (/) eRP (3rd of second dimension of summed_eRPs)
--> only those that intersect pol2 ChIP
--> --> regardless of gene group'''
plot_nogg(summed_eRPs[:,2], pol2_vals, 'sum_dc_d_gamma{}_onlyint.png'.format(gamma), np.max(summed_eRPs[:,2])+5, xlim=True)



#subselect just the gene_group, eRP value, num_selected
df_v_sub = df_v.iloc[:,[6,8,5]]
#add a signal value of 0
df_v_sub.insert(2,15,0)
#add a distance correction in the same manner as df_int_sub
df_v_sub.insert(4, 'dist_correction', pd.DataFrame.abs(((df_v.iloc[:,2]-df_v.iloc[:,1])/2) - df_v.iloc[:,3])**-gamma)
#insert ccRE_name column so that I can find the unique ccREs more easily
df_v_sub.insert(0, 'ccRE_name', df_v[0] + "_" + df_v[1].astype(str) + "_" + df_v[2].astype(str))

#combine the two
df_full = pd.concat([df_int_sub, df_v_sub])

unique_ccREs_full = np.array(df_full.ccRE_name.unique())

'''Plot everything assuming no intersection means signal value of 0 '''
summed_eRPs = np.zeros((unique_ccREs_full.shape[0], 3))
pol2_vals = np.zeros(unique_ccREs_full.shape[0])
values_gg = np.full((unique_ccREs_int.shape[0], 4, 4), np.nan) # first dimension is ccRE; second dimension is sum or pol2 val ;third dimension is gg
for i, ccRE in enumerate(unique_ccREs_full):
    '''regardless of gene group '''
    mask = df_full.iloc[:,0] == ccRE
    ''' sum of raw'''
    sum_e = pd.DataFrame.sum(df_full[mask].iloc[:,2])
    ''' sum of distance corrected '''
    sum_e_dc_p = pd.DataFrame.sum(df_full[mask].iloc[:,2]*df_full[mask].iloc[:,5])
    sum_e_dc_d = pd.DataFrame.sum(df_full[mask].iloc[:,2]/df_full[mask].iloc[:,5])
    summed_eRPs[i,[0,1,2]] = [sum_e, sum_e_dc_p, sum_e_dc_d]
    pol2_val = df_full[mask].iloc[:,3].unique()[0]
    pol2_vals[i] = pol2_val

    # '''taking into account gene group'''
    # mask1 = (df_full.iloc[:,0] == ccRE) & (df_full.iloc[:,1] == 1)
    # mask2 = (df_full.iloc[:,0] == ccRE) & (df_full.iloc[:,1] == 2)
    # mask3 = (df_full.iloc[:,0] == ccRE) & (df_full.iloc[:,1] == 3)
    # mask4 = (df_full.iloc[:,0] == ccRE) & (df_full.iloc[:,1] == 4)
    #
    # sum_e1 = pd.DataFrame.sum(df_full[mask1].iloc[:,2])
    # sum_e_dc_p1 = pd.DataFrame.sum(df_full[mask1].iloc[:,2]*df_full[mask1].iloc[:,5])
    # sum_e_dc_d1 = pd.DataFrame.sum(df_full[mask1].iloc[:,2]/df_full[mask1].iloc[:,5])
    # pol2_val = df_full[mask1].iloc[:,3].unique()[0]
    # values_gg[i, [0,1,2,3], 0] = [sum_e1, sum_e_dc_p1, sum_e_dc_d1, pol2_val]
    #
    # sum_e2 = pd.DataFrame.sum(df_full[mask2].iloc[:,2])
    # sum_e_dc_p2 = pd.DataFrame.sum(df_full[mask2].iloc[:,2]*df_full[mask2].iloc[:,5])
    # sum_e_dc_d2 = pd.DataFrame.sum(df_full[mask2].iloc[:,2]/df_full[mask2].iloc[:,5])
    # pol2_val = df_full[mask2].iloc[:,3].unique()[0]
    # values_gg[i, [0,1,2,3], 1] = [sum_e2, sum_e_dc_p2, sum_e_dc_d2, pol2_val]
    #
    # sum_e3 = pd.DataFrame.sum(df_full[mask3].iloc[:,2])
    # sum_e_dc_p3 = pd.DataFrame.sum(df_full[mask3].iloc[:,2]*df_full[mask3].iloc[:,5])
    # sum_e_dc_d3 = pd.DataFrame.sum(df_full[mask3].iloc[:,2]/df_full[mask3].iloc[:,5])
    # pol2_val = df_full[mask3].iloc[:,3].unique()[0]
    # values_gg[i, [0,1,2,3], 2] = [sum_e3, sum_e_dc_p3, sum_e_dc_d3, pol2_val]
    #
    # sum_e4 = pd.DataFrame.sum(df_full[mask4].iloc[:,2])
    # sum_e_dc_p4 = pd.DataFrame.sum(df_full[mask4].iloc[:,2]*df_full[mask4].iloc[:,5])
    # sum_e_dc_d4 = pd.DataFrame.sum(df_full[mask4].iloc[:,2]/df_full[mask4].iloc[:,5])
    # pol2_val = df_full[mask4].iloc[:,3].unique()[0]
    # values_gg[i, [0,1,2,3], 3] = [sum_e4, sum_e_dc_p4, sum_e_dc_d4, pol2_val]


#np.save('all_values_gg_gamma{}.npy'.format(gamma), values_gg)

'''Plot summed raw eRP (1st of second dimension of summed_eRPs)
--> everything assuming no intersection means signal value of 0
--> --> regardless of gene group'''
plot_nogg(summed_eRPs[:,0], pol2_vals, 'sum_raw.png', np.max(summed_eRPs[:,0])+5, xlim=True)

'''Plot summed distance corrected (*) eRP (2nd of second dimension of summed_eRPs)
--> everything assuming no intersection means signal value of 0
--> --> regardless of gene group'''
plot_nogg(summed_eRPs[:,1], pol2_vals, 'sum_dc_p_gamma{}.png'.format(gamma), np.max(summed_eRPs[:,1])+5, xlim=True)

'''Plot summed distance corrected (/) eRP (3rd of second dimension of summed_eRPs)
--> everything assuming no intersection means signal value of 0
--> --> regardless of gene group'''
plot_nogg(summed_eRPs[:,2], pol2_vals, 'sum_dc_d_gamma{}.png'.format(gamma), np.max(summed_eRPs[:,2])+5, xlim=True)
