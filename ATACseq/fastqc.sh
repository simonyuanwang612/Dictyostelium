#!/bin/bash
#SBATCH -n 1
#SBATCH -p priority
#SBATCH -t 0-00:30:00
#SBATCH -o fastqc_dir_%j.out

# invoke this script like:
# sbatch fastqc.sh sample.fastq
# where sample.fastq is an actual fastq filename 

module load fastqc/0.11.5

fastqc $1 -o .
