#!/bin/bash

bowtie -S -p 4 -a -v 0 -f --norc /data/genomes/dm3/dm3 pentamers.fasta | samtools view -Sbo pentamers.bam -
