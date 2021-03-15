#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with                                            # -N 1 means all cores will be on the same node)
#SBATCH -t 0-08:00                         # Runtime in D-HH:MM format
#SBATCH -p priority                         # Partition to run in
#SBATCH --mem=60000
module load gcc/6.2.0
module load samtools/1.3.1

samtools view -h $1 | ../../scripts/removeChrome.py  - - DDB0169550 | samtools view -b - > $2


