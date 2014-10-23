#!/bin/bash

PZ_PATH=$HOME/Dropbox/Software/peakzilla/peakzilla.py
OUTDIR=peakzilla

if [ ! -e $OUTDIR ]; then
  mkdir $OUTDIR
fi

echo Peak-calling on ChIP-exo samples:
parallel -uj 4 pypy $PZ_PATH -f 16 -l $OUTDIR/\`basename {} bed\`log {} \> $OUTDIR/\`basename {} bed\`tsv ::: *chipexo*.bed
wc -l $OUTDIR/*chipexo*.tsv

echo Peak-calling on ChIP-nexus samples:
parallel -uj 4 pypy $PZ_PATH -f 16 -l $OUTDIR/\`basename {} bed\`log {} \> $OUTDIR/\`basename {} bed\`tsv ::: *chipnexus*.bed
wc -l $OUTDIR/*chipnexus*.tsv

echo Peak-calling on ChIP-seq samples:
parallel -uj 4 pypy $PZ_PATH -l $OUTDIR/\`basename {} bed\`log {} \> $OUTDIR/\`basename {} bed\`tsv ::: *chipseq*.bed
wc -l $OUTDIR/*chipseq*.tsv

