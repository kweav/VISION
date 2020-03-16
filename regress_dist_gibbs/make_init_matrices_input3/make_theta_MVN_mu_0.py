#!/usr/bin/env python3

import numpy as np
import sys
import init_regress

#k = int(sys.argv[1]) #52

#output_file = 'theta_MVN_mu_0.npy'
output_file2 = 'init_theta_init3.npy'
#to_write = np.zeros(k)
#np.save(output_file, to_write)
#to_write = np.random.uniform(low=-10, high=10, size=(k,))
coeffs = init_regress.main()
to_write = np.hstack((coeffs, coeffs))
np.save(output_file2, to_write)
