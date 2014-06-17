#! /usr/bin/env bash

#BSUB -J rtab[1-8]
#BSUB -o rtab.%J.%I.out
#BSUB -e rtab.%J.%I.err

<<DOC
Generate 2-bin UTR windows and 20-bin CDS windows for mRNA.  Outputs tab-delimitted file with header and columns
of chrom, start, stop, count, gene, bin number (1-20), strand, and category (ex. mRNA, rRNA, tRNA)
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

results=$RESULT/$sample
bgresults=$RESULT/$sample/bedgraphs

strands=("pos" "neg")

for strand in $strands; do
    for bgfile in $(ls $bgresults/*$strand.CPMs.bg); do
        windowfile="$bgresults/$(basename $bgfile .bg).window.$strand.tab"
        fullwindowfile="$bgresults/$(basename $bgfile .bg).window.full.tab"

        # change to python script at some point in future?
        if [[ $strand == 'pos']]; then
            bedtools intersect -a $bgfile -b $MRNAWINDOWS -wao \
                | awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tcount\tgene\tbin\tstrand\tcat"} $10=="+" {split($8, a, "_"); print $1, $2, $3, $4, a[1], a[2],"+","exon"}' \
                >> $windowfile

            bedtools intersect -a $bgfile -b $UTRWINDOWS -wao \
                | awk 'BEGIN{OFS="\t"} $10=="+" {split($8,a,"_"); print $1, $2, $3, $4, a[1],a[3],"+",a[2]}' \
                >> $windowfile

        if [[ $strand == 'neg']]; then
            bedtools intersect -a $bgfile -b $MRNAWINDOWS -wao \
                | awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tcount\tgene\tbin\tstrand\tcat"} $10=="-" {split($8, a, "_"); print $1, $2, $3, $4, a[1], a[2],"-","exon"}' \
                >> $windowfile

            bedtools intersect -a $bgfile -b $UTRWINDOWS -wao \
                | awk 'BEGIN{OFS="\t"} $10=="-" {split($8,a,"_"); print $1, $2, $3, $4, a[1],a[3],"-",a[2]}' \
                >> $windowfile

        cat $windowfile >> $fullwindowfile
    done
done
