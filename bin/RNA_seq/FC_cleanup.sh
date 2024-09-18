#!/bin/bash
set -e

number=${1}
number2=${2}
counts=${3}

let sum=$number+$number2
echo $counts
echo $sum


cut -f 1,7-${sum} ${counts} > /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts2.txt 


tail -n+2 /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts2.txt > /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts2_Clean.txt
