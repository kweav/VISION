#!/bin/bash

source activate basic

date; time python simple_joint_sampler.py --where comp --threads 30 --chroms chr19 &> run_sampler.txt

conda deactivate
