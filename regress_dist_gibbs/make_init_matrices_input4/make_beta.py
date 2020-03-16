#!/usr/bin/env python3

import numpy as np
import sys
from scipy import stats

'''Usage: ./make_beta.py 52 init_theta_init4.npy init_Sigma_init4.npy'''

k = int(sys.argv[1]) #52
mean_vector = sys.argv[2]
cov_matrix = sys.argv[3]

init_beta = stats.multivariate_normal.rvs(mean=np.load(mean_vector), cov=np.load(cov_matrix))
print(init_beta.shape[0] == k)
to_write = 'init_beta_init4.npy'
np.save(to_write, init_beta)
