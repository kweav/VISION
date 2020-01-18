#!/usr/bin/env python3

import numpy as np
import sys

'''Usage: ./any_all_0.py TSS_expression.npz '''

exp_file = sys.argv[1]


npzfile = np.load(exp_file, allow_pickle = True)
exp_values = npzfile['exp']
print(np.sum(np.sum(exp_values, axis=1) == 0))
