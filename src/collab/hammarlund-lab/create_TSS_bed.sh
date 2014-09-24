#! /usr/bin/env bash

<<DOC
Grab TSS from ce10 refGene BED12 format, output TSSs as single bp in BED6 format
DOC

genomedir="/vol2/home/speach/ref/genomes/ce10/"
refbed="${genomedir}refGene.ce10.bed"
tmptssbed="${genomedir}refGene.ce10.TSS.unsorted.bed"
tssbed="${genomedir}refGene.ce10.TSS.bed"

# TSS of positive strand genes
awk 'BEGIN{OFS="\t"} $6 == "+" {print $1, $7, $7+1, $4, $5, $6}' $refbed \
    > $tmptssbed

# TSS of negative strand genes
awk 'BEGIN{OFS="\t"} $6 == "-" {print $1, $7, $7+1, $4, $5, $6}' $refbed \
    >> $tmptssbed

# Sort bed, remove temporary bed
sort -k1,1 -k2,2n $tmptssbed > $tssbed
rm -rf $tmptssbed