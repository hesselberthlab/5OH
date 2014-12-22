#! /usr/bin/env bash

<<DOC
Generate user-friendly summary of changepoint hits
DOC

for file in $( ls *changepoints.bed ); do
    sample=$(basename $file .changepoints.bed)
    awk 'BEGIN{OFS="\t"} {print $4, $5, $1:$2}' $file \
    | sort -k2,2gr \
    > ${sample}.results.tab
done