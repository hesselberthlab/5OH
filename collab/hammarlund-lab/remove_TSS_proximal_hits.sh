#! /usr/bin/env bash

<<DOC
Calculate distance from changepoint to nearest TSS.  Remove hits within a given distance.

Usage: remove_TSS_proximal_hits.sh $compbed $tssdistance
DOC

#Read in arguments
compbed=$1
tssdistance=$2

compname=$(basename $compbed .bed)
output=${compname}.${tssdistance}-tss.bed

genomedir="/vol2/home/speach/ref/genomes/ce10/"
tssbed="${genomedir}refGene.ce10.TSS.bed"

closestBed -a $compbed -b $tssbed -s \
    | awk -v tssd=$tssdistance 'BEGIN{OFS="\t"} ($2-$8) > tssd || ($8-$2) > tssd {print $1, $2, $3, $4, $5, $6, $2-$8}' \
    > $output