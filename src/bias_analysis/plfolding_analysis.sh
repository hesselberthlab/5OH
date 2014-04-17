#! /usr/bin/env bash

# example/skeleton of plfolding_analysis pipeline
# i'll probably just write everything in python all the time, eventually

# Run once only
# python ~/devel/50H/src/plfolding_analysis.py sacCer1MRNA_s.fa > sacCer1MRNA_bpProb.bed

samplebg=~/projects/5OH/results/20140107/*S4.pos.bg
probbed=sacCer1MRNA_bpProb.bed

bedtools intersect -a $samplebg -b $probbed -wao \
	> pl_intersect_raw.tab

awk 'BEGIN{print "count\tprob\tgene"; OFS="\t"} $11 == 1 && $10 == "+" {print $4, $8, $9}' pl_intersect_raw.tab \
	> pl_intersect.tab

Rscript plfolding_analysis.R
