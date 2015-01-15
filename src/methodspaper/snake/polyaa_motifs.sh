#! /usr/bin/env bash

<<DOC
Generate csv files for R analysis of polyaa motifs and make meta plots via
R script.

Usage (Snakemake):
rule motifs:
    input: RESULT + "{sample}." + ASAL_TEXT + ".bam"
    params: POLYALL, MOTIFS, POLYAACODEDIR
    output: MOTIFS + "{sample}." + ASAL_TEXT + ".poly-all.avgs.csv", MOTIFS + "{sample}." + ASAL_TEXT + ".poly-all.pdf"
    shell: "bash polyaa_motifs.sh {input} {params} {output}"
DOC

set -o nounset -o pipefail -o errexit -x

# Read variables from command line
bam=$1
polyall=$2
outputdir=$3
codedir=$4
csvfile=$5
outputpdf=$6

samplename=$(basename $bam .bam)
sample=${outputdir}${samplename}
allbaseposbg=${sample}.allbase.pos.bg

# Currently interested in looking at all motifs,
# but legacy code has this as a variable
aa="-all"

# if allbase.bgs have not been created for sample, make it happen.
if [[ ! -f $allbaseposbg ]]; then
    cp $bam $3/.
    bash $codedir/create_allbase_bg.sh $sample    
fi

intersectfile=${sample}.poly${aa}.intersect.tab
CovPMpos=${sample}.allbase.CovPM.pos.bg
CovPMneg=$sample.allbase.CovPM.neg.bg

# Intersect poly-aa with allbasebg and generate csvfile for input to R script
if [[ ! -f $intersectfile ]]; then
    awk '$6 == "+"' $polyall \
        | bedtools map -a - -b $CovPMpos -c 4 -o collapse \
        > $intersectfile

    awk '$6 == "-"' $polyall \
        | bedtools map -a - -b $CovPMneg -c 4 -o collapse \
        >> $intersectfile

    python ${codedir}/polyaa_matrix.py $intersectfile $csvfile > $csvfile.motifs
fi

# Run R script to output metadata tables of all motifs
Rscript ${codedir}/plot_poly_metadata_HT.R $csvfile $outputpdf
