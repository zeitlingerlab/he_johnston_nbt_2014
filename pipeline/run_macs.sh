#!/bin/bash

SAMPLE_NAME=`basename $1 .bed`

macs2 callpeak -g dm --keep-dup=all --outdir macs -n $SAMPLE_NAME --call-summits -f BED -t $1
