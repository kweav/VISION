#!/bin/bash

source activate basic
date; time python use_clustered_regress.py --where comp --chroms $1
conda deactivate
