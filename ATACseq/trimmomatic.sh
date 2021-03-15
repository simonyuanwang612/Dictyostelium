#!/bin/bash
#SBATCH -n 4
#SBATCH -o trimmomatic_%J.out
#SBATCH -p priority
#SBATCH --mem=12G
#SBATCH -t 0-5
#SBATCH -e trimmomatic_%j.err                 
module load trimmomatic/0.36

java -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE \
-threads 4 -phred33 \
../../data/"$1"_1.fq.gz ../../data/"$1"_2.fq.gz \
"$1".R1.trimmed.paired.fastq.gz "$1".R1.trimmed.unpaired.fastq.gz \
"$1".R2.trimmed.paired.fastq.gz "$1".R2.trimmed.unpaired.fastq.gz \
ILLUMINACLIP:$TRIMMOMATIC/adapters/NexteraPE-PE.fa:2:30:10:2:true  \
LEADING:8 \
TRAILING:8 \
MINLEN:30
