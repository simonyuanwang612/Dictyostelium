#!/bin/bash

#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with                                            # -N 1 means all cores will be on the same node)
#SBATCH -t 0-08:00                         # Runtime in D-HH:MM format
#SBATCH -p short                         # Partition to run in
#SBATCH --mem=60000
module load gcc/6.2.0
module load python/2.7.12
module load macs2/2.1.1.20160309

macs2 callpeak -t $1 -f BAMPE -n $2 -g 3.4e7 --keep-dup all --cutoff-analysis --outdir $3 
