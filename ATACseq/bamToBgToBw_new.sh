#!/bin/bash                         

file=$1
echo $file
genome=$2
echo $genome

script=/ye/netapp/jimmie.ye/tools/lab.scripts/bamToBgToBw

#cd ${file}/bams
# bam to bedgrah
bedtools genomecov -ibam ${file}.unique.bam -bg | sort -k1,1 -k2,2n > ${file}.bg

# bedgraph to bigwig
${script}/bedGraphToBigWig ${file}.bg ${script}/${genome}.chrom.sizes ${file}.bw
