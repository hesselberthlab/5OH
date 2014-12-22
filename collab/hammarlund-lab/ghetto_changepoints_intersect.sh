#! /usr/bin/env bash

<<DOC
Intersect stranded bedgraphs with ensemble genes, calculate change points, and
output BED formatted file where column 5 is score ratio.

todo: make highthruput and fix stuff
DOC

pos_intersect=RTCB_TM_012_CTTGTA_L002.strand.pos.intersect.tab
neg_intersect=RTCB_TM_012_CTTGTA_L002.strand.neg.intersect.tab
output=RTCB_TM_012_CCTGTA_L002.changepoints.tab

# positive strand 
awk 'BEGIN{OFS="\t"} $6 == "+"' ensGene.ce10.bed \
    | cut -f 1-6 \
    | bedtools intersect -a - -b RTCB_TM_012_CTTGTA_L002.strand.pos.bg -wb \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $10, $6}' \
    > $pos_intersect

# negative strand
awk 'BEGIN{OFS="\t"} $6 == "-"' ensGene.ce10.bed \
    | cut -f 1-6 \
    | bedtools intersect -a - -b RTCB_TM_012_CTTGTA_L002.strand.neg.bg -wb \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $10, $6}' \
    > $neg_intersect

python ghetto_changepoints.py -p $pos_intersect -n $neg_intersect > $output



