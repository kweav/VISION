#!/bin/bash

source activate basic

date; time python simple_joint_sampler.py --where comp --threads 30 --chroms chr1 &> run_simple_sampler_chr1.txt

conda deactivate
