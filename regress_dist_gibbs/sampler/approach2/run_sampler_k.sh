#!/bin/bash

source activate basic

date; time python joint_sampler_better_k.py --where comp --threads 30 --chroms chr1 --init_k 0.8 > run_sampler_k_chr1.txt

conda deactivate
