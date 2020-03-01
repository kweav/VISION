#!/bin/bash

beta_a=5
beta_b=2
num_points=50
step=0.01
x_lim0=0.01
x_lim1=1

date;time ./approx_beta.py $beta_a $beta_b $num_points $x_lim0 $x_lim1 $step > outfit.txt
