# Kallisto indexing

kallisto index dicty_chromosomal.fasta --index=dicty_base_idx

# Kallisto quanting
sbatch -p short -n 4 --mem 5G -t 0-2 -o lib13.base.log wrap="kallisto quant -i dicty_base_idx -o LIB13.base.out --bootstrap-samples=100 ../../../data/LIB13_R1.fastq.bz2  ../../../data/LIB13_R2.fastq.bz" 
