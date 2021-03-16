#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with                                            # -N 1 means all cores will be on the same node)
#SBATCH -t 0-08:00                         # Runtime in D-HH:MM format
#SBATCH -p priority                          # Partition to run in
#SBATCH --mem=100000
module load gcc/6.2.0
module load samtools/1.3.1
module load bowtie2/2.2.9

bowtie2 -x dicty_idx -1 $1 -2 $2 | samtools view -bS - | samtools sort - > $3


