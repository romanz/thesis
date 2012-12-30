#!/usr/bin/env python

import sys

script = '''#!/bin/sh
#PBS -N ElectroKinetics
#PBS -q all_l_p
#PBS -M  romanz@cs.technion.ac.il
#PBS -mbea
#PBS -l select=1:ncpus=4

PBS_O_WORKDIR=$HOME/Code/thesis/src/
cd $PBS_O_WORKDIR

echo "%s" | matlab -nodisplay -nosplash 
'''
args = sys.argv[1]
print script % args
