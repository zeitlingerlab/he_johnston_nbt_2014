#!/bin/bash

# 1st parameter: Preprocessed FASTQ
# 2nd parameter: Path to bowtie reference genome

bowtie -S -p 3 --chunkmbs 512 -k 1 -m 1 -v 2 --best --strata \
       $2 <(cutadapt -m 22 -O 4 -e 0.2 \
       -a AGATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
       -a GATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
       -a ATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
       -a TCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
       -a CGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
       $1) | samtools view -Sb -F 4 -o `basename $1 _processed.fastq.gz`.bam -
 
 
