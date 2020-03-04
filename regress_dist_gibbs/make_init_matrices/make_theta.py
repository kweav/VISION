#!/usr/bin/env python3

'''Note this matrix must be positve definite
    check symmetric and cholesky decomposition'''

import numpy as np
import sys
from scipy import linalg

def is_pos_def(x):
    chol = linalg.cholesky(x)
    #will raise error if can't decompose
    return x.T == x


output_file = 'Sigma_invwishart_S_0.npy'
output_file2 = 'init_Sigma.npy'

print(is_pos_def(to_write), flush=True)
