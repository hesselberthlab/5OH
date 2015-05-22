#! /usr/bin/env bash

#BSUB -J windows[1-8]
#BSUB -o windows.%J.%I.out
#BSUB -e windows.%J.%I.err

<<DOC
Generate 2-bin UTR windows and 20-bin CDS windows for mRNA.  Outputs tab-delimitted file with header and columns
of chrom, start, stop, count, gene, bin number (1-20), strand, and category (ex. mRNA, rRNA, tRNA)
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

for bgfile in $(ls $bgresults/*both.CPMs.bg); do
    # Parse basename of file and create intersection files
    posbg="$bgresults/$(basename $bgfile .both.CPMs.bg).pos.CPMs.bg"
    negbg="$bgresults/$(basename $bgfile .both.CPMs.bg).neg.CPMs.bg"
    
    rfile="$rtables/$(basename $bgfile .bg).windows.tab"

    # change to python script at some point in future?
    bedtools intersect -a $posbg -b $MRNAWINDOWS -wao \
        | awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tcount\tgene\tbin\tstrand\tcat"} $10=="+" {split($8, a, "_"); print $1, $2, $3, $4, a[1], a[2],"+","exon"}' \
        > $rfile

    bedtools intersect -a $posbg -b $UTRWINDOWS -wao \
        | awk 'BEGIN{OFS="\t"} $10=="+" {split($8,a,"_"); print $1, $2, $3, $4, a[1],a[3],"+",a[2]}' \
        >> $rfile

    bedtools intersect -a $negbg -b $MRNAWINDOWS -wao \
        | awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tcount\tgene\tbin\tstrand\tcat"} $10=="-" {split($8, a, "_"); print $1, $2, $3, $4, a[1], a[2],"-","exon"}' \
        >> $rfile

    bedtools intersect -a $negbg -b $UTRWINDOWS -wao \
        | awk 'BEGIN{OFS="\t"} $10=="-" {split($8,a,"_"); print $1, $2, $3, $4, a[1],a[3],"-",a[2]}' \
        >> $rfile
done
