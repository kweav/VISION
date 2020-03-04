#!/usr/bin/env python3

import numpy as np
import sys

k = int(sys.argv[1]) #52

output_file = 'theta_MVN_mu_0.npy'
output_file2 = 'init_theta.npy'
to_write = np.zeros(k)
np.save(output_file, to_write)
np.save(output_file2, to_write)
