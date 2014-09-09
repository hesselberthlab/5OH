#! /usr/bin/env bash

<<DOC
Parse genes into 20 mRNA bins and 2 UTR bins.

input: RESULT + "{sample}.align." + ALIGN_MODE + ".strand.pos.CPMs.bg",
       RESULT + "{sample}.align." + ALIGN_MODE + ".strand.neg.CPMs.bg"
params: MRNAWINDOWS, UTRWINDOWS
output: GRAPHS + "{sample}.align." + ALIGN_MODE + ".windows.tab"
DOC

set -o nounset -o pipefail -o errexit -x

# Read variables from command line
# TODO: MRNAWINDOWS and UTRWINDOWS are currently read in as variables
# These can be dynamically created in this code, allowing using to
# specify bins or UTR inclusion in Snakefile configuration
posbg=$1
negbg=$2
MRNAWINDOWS=$3
UTRWINDOWS=$4
rfile=$5


# Intersect bedgraph files with MRNA and UTR window files to generate R-compatible table
# TODO: (Consider changing to Python at some point to make other people happy.)
bedtools intersect -a $posbg -b $MRNAWINDOWS -wao \
    | awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tcount\tgene\tbin\tstrand\tcat"} $10=="+" {split($8, a, "_"); print $1, $2, $3, $4, a[1], a[2],"+","exon"}' \
    > $rfile

bedtools intersect -a $posbg -b $UTRWINDOWS -wao \
    | awk 'BEGIN{OFS="\t"} $10=="+" {split($8,a,"_"); print $1, $2, $3, $4, a[1],a[3],"+",a[2]}' \
    >> $rfile

bedtools intersect -a $negbg -b $MRNAWINDOWS -wao \
    | awk 'BEGIN{OFS="\t"} $10=="-" {split($8, a, "_"); print $1, $2, $3, $4, a[1], a[2],"-","exon"}' \
    >> $rfile

bedtools intersect -a $negbg -b $UTRWINDOWS -wao \
    | awk 'BEGIN{OFS="\t"} $10=="-" {split($8,a,"_"); print $1, $2, $3, $4, a[1],a[3],"-",a[2]}' \
    >> $rfile
