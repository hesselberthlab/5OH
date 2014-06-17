#! /usr/bin/env bash

#BSUB -J rtab[1-8]
#BSUB -o rtab.%J.%I.out
#BSUB -e rtab.%J.%I.err

<<DOC
Generate tables for R analysis and graphing
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

results=$RESULT/$sample
bgresults=$RESULT/$sample/bedgraphs

strands=("pos" "neg")

for strand in $strands; do
    for bgfile in $(ls $bgresults/*$strand.CPMs.bg); do
        rfile="$bgresults/$(basename $bgfile .bg).intersect.tab"

        if [[ $strand == 'pos']]; then
              bedtools intersect -a $bgfile -b $FULLGFF -wao \
              | awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tcount\tcat\tgene\tsite"} $10=="+" {split($8, a, ":"); print $1, $2, $3, $4, a[1], a[2], $1 ":" $2}' \
              > $rfile

        if [[ $strand == 'neg']]; then
              bedtools intersect -a $bgfile -b $FULLGFF -wao \
              | awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tcount\tcat\tgene\tsite"} $10=="-" {split($8, a, ":"); print $1, $2, $3, $4, a[1], a[2], $1 ":" $2}' \
              > $rfile
    done
done
