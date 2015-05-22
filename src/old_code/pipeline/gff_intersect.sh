#! /usr/bin/env bash

<<DOC
Generate tables for R analysis and graphing

input: RESULT + "{sample}.align." + ALIGN_MODE + ".strand.pos.CPMs.bg",
       RESULT + "{sample}.align." + ALIGN_MODE + ".strand.neg.CPMs.bg"
params: FULLGFF
output: RESULT + "{sample}.align." + ALIGN_MODE + ".gffintersect.tab"
DOC

set -o nounset -o pipefail -o errexit -x

# Read variables from command line
posbg=$1
negbg=$2
FULLGFF=$3
rfile=$4

# intersect bg file with master GFF; parse appropriately
touch $rfile

bedtools intersect -a $posbg -b $FULLGFF -wao \
    | awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tcount\tcat\tgene\tsite"} $10=="+" {split($8, a, ":"); print $1, $2, $3, $4, a[1], a[2], $1 ":" $2}' \
    > $rfile

bedtools intersect -a $negbg -b $FULLGFF -wao \
    | awk 'BEGIN{OFS="\t"} $10=="-" {split($8, a, ":"); print $1, $2, $3, $4, a[1], a[2], $1 ":" $2}' \
    >> $rfile
