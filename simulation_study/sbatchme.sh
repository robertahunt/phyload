#!/bin/bash
#SBATCH
#SBATCH -o sbatchme.out
#SBATCH -e sbatchme.err
##SBATCH -p largenode
##SBATCH -c 24
##SBATCH --mem=30000

revpath=/home/wdewitt/wdewitt_fast/revbayes/projects/cmake/rb
source activate phyload
scons --revpath=${revpath} -j 70 -k
