#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with                                            # -N 1 means all cores will be on the same node)
#SBATCH -t 0-08:00                         # Runtime in D-HH:MM format
#SBATCH -p priority                          # Partition to run in
#SBATCH --mem=40000                          # Memory total in MB (for all cores)
module load java/jdk-1.8u112 
module load picard/2.8.0 

java -jar $PICARD/picard-2.8.0.jar MarkDuplicates INPUT=$1 OUTPUT=$2 METRICS_FILE=$3 ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT



