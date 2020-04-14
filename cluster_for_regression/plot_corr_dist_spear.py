#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

npz_file = np.load(sys.argv[1])
saved_matrix = npz_file['corrs_and_dists'] #3 dim: 1 TSS, 2 CRE, 3rd dim 0 index is little d distance, 1 index is Cap D distance assuming only two clustered, 2 index is spearmanr
fig, (ax1, ax2) = plt.subplots(1,2)
ax1.set_title('{} d_ij vs spearmanr_ij'.format(sys.argv[2]))
ax1.set_xlabel('spearmanr_ij')
ax1.set_ylabel('d_ij')
ax2.set_title('{} D_ij vs spearmanr_ij'.format(sys.argv[2]))
ax2.set_xlabel('spearmanr_ij')
ax2.set_ylabel('D_ij')
isnan_bool = ~np.isnan(saved_matrix[:,:,2].reshape(-1,))
print(sys.argv[2], isnan_bool.shape[0]-np.sum(isnan_bool), flush=True)
ax1.scatter(saved_matrix[:,:,2].reshape(-1,)[isnan_bool], saved_matrix[:,:,0].reshape(-1,)[isnan_bool])
pearson1, pval = stats.pearsonr(saved_matrix[:,:,2].reshape(-1,)[isnan_bool], saved_matrix[:,:,0].reshape(-1,)[isnan_bool])
ax1.annotate("Pearson R: {0:.3f}".format(pearson1), (1,1))
print("Pearson1", sys.argv[2], pearson1, flush=True)
ax2.scatter(saved_matrix[:,:,2].reshape(-1,)[isnan_bool], saved_matrix[:,:,1].reshape(-1,)[isnan_bool])
pearson2, pval = stats.pearsonr(saved_matrix[:,:,2].reshape(-1,)[isnan_bool], saved_matrix[:,:,1].reshape(-1,)[isnan_bool])
ax2.annotate("Pearson R: {0:.3f}".format(pearson2), (1,1))
print("Pearson2", sys.argv[2], pearson2, flush=True)
fig.savefig('{}_dist_vs_spear.png'.format(sys.argv[2]))
plt.close(fig)
