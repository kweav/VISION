#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
input_file, output_file = sys.argv[1:]

df = pd.read_csv(input_file, delimiter='\t')
df_gene_len = df.iloc[:,2] - df.iloc[:,1]
df_counts = df.iloc[:,6:]
df_summed = df_counts.sum(axis=0).transpose()

df_tpm = df_counts.div(df_gene_len, axis=0).div(df_summed, axis=1).mul(1000000)
df_log2_tpm = np.log2(df_tpm + 0.001)
df_out = pd.concat((df.iloc[:,0:6], df_log2_tpm), axis=1)
df_out.to_csv(output_file, sep='\t', index=False, float_format='%.3f')
