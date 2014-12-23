#! /usr/bin/env bash

<<DOC
Creating the polyaa matrix CSV for useage in R script for polyaa meta data
heatmap plots
DOC

codedir=~/devel/5OH/src/codon_analysis
resultsdir=/vol2/home/speach/projects/5OH/results/methodspaper/human
aa_arg="-all"
aa="all"
polybed=${resultsdir}/hg19.poly${aa}_genes.bed

#python ${codedir}/get_poly_stretches_location_mismatch.py $aa_arg \
#    | sort -k1,1 -k2,2n \
#    > $polybed

sample=SP30.assembly-canonicalhg19.align-uniq
sampleprefix=${resultsdir}/${sample}
bam=${sampleprefix}.bam
allbasebg=${sampleprefix}.allbase.bg

# for the human stuff, only positive strand matters
bedtools genomecov -ibam $bam -d -strand + \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $2 + 1, $3}' \
    > $allbasebg

intersectfile=${sampleprefix}.poly${aa}.intersect.tab
csvfile=${sampleprefix}.poly${aa}.avgs.csv

bedtools map -a $polybed -b $allbasebg -c 4 -o collapse \
    > $intersectfile

python ${codedir}/polyaa_matrix.py $intersectfile $csvfile \
    > $csvfile.motifs
