#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
input_file, output_file = sys.argv[1:]

'''TPM is defined as the following:
A* 1/sum(A)*10^6 where A is the reads per kilobase or A = total reads mapped to gene/gene length in kb
Notice the order of operations
'''

df = pd.read_csv(input_file, delimiter='\t')
df_gene_len = (df.iloc[:,2] - df.iloc[:,1])/1000 #finding gene length in kb
df_counts = df.iloc[:,6:]
df_rpk = df_counts.div(df_gene_len, axis=0) #reads per kilobase
df_summed_rpk = df_rpk.sum(axis=0)

df_tpm = df_rpk.div(df_summed_rpk, axis=1).mul(1000000)
df_log2_tpm = np.log2(df_tpm + 0.001)
df_out = pd.concat((df.iloc[:,0:6], df_log2_tpm), axis=1)
df_out.to_csv(output_file, sep='\t', index=False, float_format='%.3f')
