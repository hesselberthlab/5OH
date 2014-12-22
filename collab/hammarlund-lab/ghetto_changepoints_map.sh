#! /usr/bin/env bash

<<DOC
Map stuff and calc change points.
DOC

pos_map=RTCB_TM_012_CTTGTA_L002.strand.pos.mapping.tab
neg_map=RTCB_TM_012_CTTGTA_L002.strand.neg.mapping.tab

# positive strand 
awk 'BEGIN{OFS="\t"} $6 == "+"' ensGene.ce10.bed \
    | cut -f 1-6 \
    | bedtools map -a - -b RTCB_TM_012_CTTGTA_L002.strand.pos.bg -c 4 -o collapse \
    > $pos_map

# negative strand
awk 'BEGIN{OFS="\t"} $6 == "-' ensGene.ce10.bed \
    | cut -f 1-6 \
    | bedtools map -a - -b RTCB_TM_012_CTTGTA_L002.strand.neg.bg -c 4 -o collapse \
    > $neg_map

python ghetto_changepoints.py $pos_map $neg_map



