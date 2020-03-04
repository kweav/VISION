#!/usr/bin/env python3

'''Note this matrix must be positve definite
    check symmetric and cholesky decomposition'''

import numpy as np
import sys
from scipy import linalg

def is_pos_def(x):
    chol = linalg.cholesky(x)
    #will raise error if can't decompose
    return np.allclose(x, x.T)

k = int(sys.argv[1]) #52

sign_dict = {1: {2: -1, #note that the numbers in the sign_dict are based on actual IDEAS state not loc. Recall no 0, so everything is shifted IDEAS_state - 1
                 3: -1,
                 8: +1,
                 14: +1,
                 16: -1,
                 17: +1,
                 20: -1,
                 22: -1,
                 25: +1},
             2: {1: -1,
                 3: +1,
                 5: -1,
                 6: -1,
                 8: -1,
                 10: -1,
                 12: -1,
                 14: -1,
                 16: +1,
                 17: -1,
                 18: -1,
                 19: -1,
                 20: +1,
                 21: -1,
                 22: +1,
                 23: -1,
                 24: -1,
                 25: -1,
                 26: -1},
            3: {1: -1,
                2: +1,
                5: -1,
                6: -1,
                8: -1,
                10: -1,
                12: -1,
                14: -1,
                16: +1,
                17: -1,
                18: -1,
                19: -1,
                20: +1,
                21: -1,
                22: +1,
                23: -1,
                24: -1,
                25: -1,
                26: -1},
            4:{5: +1},
            5:{2: -1,
               3: -1,
               4: +1,
               16: -1,
               20: -1,
               22: -1},
            6:{2: -1,
               3: -1,
               10: +1,
               11: +1,
               12: +1,
               15: +1,
               16: -1,
               18: +1,
               19: +1,
               20: -1,
               21: +1,
               22: -1,
               23: +1,
               24: +1},
            7:{9: +1,
               13: +1,
               26: +1},
            8:{1: +1,
               2: -1,
               3: -1,
               14: +1,
               16: -1,
               17: +1,
               20: -1,
               22: -1,
               25: +1},
            9:{7: +1,
               13: +1,
               26: +1},
            10:{2: -1,
               3: -1,
               6: +1,
               11: +1,
               12: +1,
               15: +1,
               16: -1,
               18: +1,
               19: +1,
               20: -1,
               21: +1,
               22: -1,
               23: +1,
               24: +1},
            11:{6: +1,
                10: +1,
                12: +1,
                15: +1,
                18: +1,
                19: +1,
                21: +1,
                23: +1,
                24: +1},
            12:{2: -1,
               3: -1,
               6: +1,
               10: +1,
               11: +1,
               15: +1,
               16: -1,
               18: +1,
               19: +1,
               20: -1,
               21: +1,
               22: -1,
               23: +1,
               24: +1},
            13:{7: +1,
                9: +1,
               26: +1},
            14:{1: +1,
                2: -1,
                3: -1,
                8: +1,
                16: -1,
                17: +1,
                20: -1,
                22: -1,
                25: +1},
            15:{6: +1,
                10: +1,
                11: +1,
                12: +1,
                18: +1,
                19: +1,
                21: +1,
                23: +1,
                24: +1},
            16:{1: -1,
                2: +1,
                3: +1,
                5: -1,
                6: -1,
                8: -1,
                10: -1,
                12: -1,
                14: -1,
                17: -1,
                18: -1,
                19: -1,
                20: +1,
                21: -1,
                22: +1,
                23: -1,
                24: -1,
                25: -1,
                26: -1},
            17:{1: +1,
                2: -1,
                3: -1,
                8: +1,
                14: +1,
                16: -1,
                20: -1,
                22: -1,
                25: +1},
            18:{2: -1,
               3: -1,
               6: +1,
               10: +1,
               11: +1,
               12: +1,
               15: +1,
               16: -1,
               19: +1,
               20: -1,
               21: +1,
               22: -1,
               23: +1,
               24: +1},
            19:{2: -1,
               3: -1,
               6: +1,
               10: +1,
               11: +1,
               12: +1,
               15: +1,
               16: -1,
               18: +1,
               20: -1,
               21: +1,
               22: -1,
               23: +1,
               24: +1},
            20:{1: -1,
                2: +1,
                3: +1,
                5: -1,
                6: -1,
                8: -1,
                10: -1,
                12: -1,
                14: -1,
                16: +1,
                17: -1,
                18: -1,
                19: -1,
                21: -1,
                22: +1,
                23: -1,
                24: -1,
                25: -1,
                26: -1},
            21:{2: -1,
               3: -1,
               6: +1,
               10: +1,
               11: +1,
               12: +1,
               15: +1,
               16: -1,
               18: +1,
               19: +1,
               20: -1,
               22: -1,
               23: +1,
               24: +1},
            22:{1: -1,
                2: +1,
                3: +1,
                5: -1,
                6: -1,
                8: -1,
                10: -1,
                12: -1,
                14: -1,
                16: +1,
                17: -1,
                18: -1,
                19: -1,
                20: +1,
                21: -1,
                23: -1,
                24: -1,
                25: -1,
                26: -1},
            23:{2: -1,
               3: -1,
               6: +1,
               10: +1,
               11: +1,
               12: +1,
               15: +1,
               16: -1,
               18: +1,
               19: +1,
               20: -1,
               21: +1,
               22: -1,
               24: +1},
            24:{2: -1,
               3: -1,
               6: +1,
               10: +1,
               11: +1,
               12: +1,
               15: +1,
               16: -1,
               18: +1,
               19: +1,
               20: -1,
               21: +1,
               22: -1,
               23: +1},
            25:{1: +1,
                2: -1,
                3: -1,
                8: +1,
                14: +1,
                16: -1,
                17: +1,
                20: -1,
                22: -1},
            26:{2: -1,
                3: -1,
                7: +1,
                9: +1,
               13: +1,
               16: -1,
               20: -1,
               22: -1}}

output_file = 'Sigma_invwishart_S_0.npy'
output_file2 = 'init_Sigma.npy'
to_write = np.zeros((k, k))
np.fill_diagonal(to_write, 25)
minimal_cov = float(sys.argv[2]) #0.25
for key in sign_dict:
    row_val = key-1
    row_val2 = int((key-1)+(k/2))
    for subkey in sign_dict[key]:
        column_val = subkey-1
        column_val2 = int((subkey-1)+(k/2))
        to_write[row_val, column_val] = sign_dict[key][subkey]*minimal_cov
        to_write[row_val2, column_val2] = sign_dict[key][subkey]*minimal_cov

np.save(output_file, to_write)
np.save(output_file2, to_write)

print(is_pos_def(to_write), flush=True)
