#!/usr/bin/env python3

'''Note this matrix must be positive semi definite
clearly symmetric so only checking for non-negative eigvals '''

import numpy as np
import sys
from scipy import linalg

def is_pos_sem_def(x):
    return np.all(linalg.eigvals(x) >= 0)

k = int(sys.argv[1]) #52

output_file = 'theta_MVN_Lambda_0.npy'
to_write = np.zeros((k,k))
np.fill_diagonal(to_write, 25)
np.save(output_file, to_write)

print(is_pos_sem_def(to_write), flush=True)
