#! /usr/bin/env bash

<<DOC
Compare data to PARS scores from Kertesz/Segal, Nature 2010 RNA structure
mapping paper.

Output table:
chrom, start, stop, CPMs, PARS score

Useage (Snakemake):
rule compare_pars:
    input: RESULT + "{sample}." + ASAL_TEXT + ".strand.all.CPMs.bg"
    params: PARSBG, RSCRIPTS
    output: PARSDIR + "{sample}." + ASAL_TEXT + ".strand.all.CPMs-PARS.cutoff-1.pdf"
    shell: "bash compare_pars.sh {input} {params} {output}
DOC

set -o nounset -o pipefail -o errexit -x

# Read variables from command line
samplebg=$1
parsbg=$2
rscripts=$3
intersecttab=$4
#intersecttab=$(dirname $pdf)$(basename $samplebg .bg).tab

bedtools intersect -a $samplebg -b $parsbg -wao \
    | awk 'BEGIN{OFS="\t"; print "chrom\tstart\tstop\tCPMs\tPARS"} $8 != "." {print $1, $2, $3, $4, $8}' \
    > $intersecttab

Rscript ${rscripts}graph_PARS_vs_CPMs.R $intersecttab
