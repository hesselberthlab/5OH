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
rtables=$RESULT/$sample/rtables
if [[ ! -d $rtables ]]; then
    mkdir -p $rtables
fi

strands=("pos" "neg")

for bgfile in $(ls $bgresults/*both.CPMs.bg); do
    # Parse basename of file and create intersection files
    posbg="$bgresults/$(basename $bgfile .both.CPMs.bg).pos.CPMs.bg"
    negbg="$bgresults/$(basename $bgfile .both.CPMs.bg).neg.CPMs.bg"
    
    rfile="$rtables/$(basename $bgfile .bg).intersect.tab"

    # intersect bg file with master GFF; parse appropriately
    bedtools intersect -a $posbg -b $FULLGFF -wao \
        | awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tcount\tcat\tgene\tsite"} $10=="+" {split($8, a, ":"); print $1, $2, $3, $4, a[1], a[2], $1 ":" $2}' \
        > $rfile

    bedtools intersect -a $negbg -b $FULLGFF -wao \
        | awk 'BEGIN{OFS="\t"} $10=="-" {split($8, a, ":"); print $1, $2, $3, $4, a[1], a[2], $1 ":" $2}' \
        >> $rfile
done
