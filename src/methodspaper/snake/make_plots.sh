#! /usr/bin/env bash

<<DOC
WORK IN PROGRESS
Run several Rscripts to generate plots of 5OH data

input: RGRAPHS + "{sample}.align." + ALIGN_MODE + ".gffintersect.tab",
       RGRAPHS + "{sample}.align." + ALIGN_MODE + ".windows.tab"
params: RSCRIPTS, RGRAPHS
output: RGRAPHS + "{sample}.align." + ALIGN_MODE + ".placeholder.txt" (output TBD!)
DOC

set -o nounset -o pipefail -o errexit -x

# Read variables from command line
intersecttab=$1
windowstab=$2
RSCRIPTS=$3
RGRAPHS=$4
output=$5
fullsample=$(basename $intersecttab .gffintersect.tab)      # Full sample with alignment info
sample=$(basename $fullsample .align.*)                     # Simple sample name. TODO: parse systematic sample name to useful one
sampledescr="temporary"                                     # TODO: pass sample descr to this script

################################
# DIFFUSE/DISCRETE SCATTER PLOTS
# These plots look subpar at the moment.
#
# minimum number of counts
mincounts=(5 10 25 50 100)

for count in "${mincounts[@]}"; do
    Rscript "${RSCRIPTS}/diffuse_discrete_scatter.R" $intersecttab $sample $sampledescr $count $RGRAPHS
done
###############################

###############################
# mRNA GENE BINNING
#
# minimum number of counts
mincounts=(5 10 25 50 100)

for count in "${mincounts[@]}"; do
    Rscript "${RSCRIPTS}/window_analysis_utrs_and_exons.R" $windowstab $sample $count $RGRAPHS
done
###############################

###############################
# PIE CHARTS OF SPECIES CAPTURED?
#
# potentially useful for paper
#
# [Rscript ...]
#
###############################


###############################
# FEATURES ...
#
# Make compatible with Jaycode
#
# [Rscript ...]
#
###############################

echo "highfive" > $output
