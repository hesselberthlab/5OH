#! /usr/bin/env bash

# Testing command line usage of R scripts

SAMPLES=(GCCTAA_S1
        TGGTCA_S2
        CACTGT_S3
        ATTGGC_S4)
DESCR=("del-XRN1,+Tm"
        "del-XRN1,DMSO"
        "WT,+Tm"
        "WT,DMSO")

MINCOUNTS=(5 10 25 50 100)

webpath=$HOME/public_html/projects/5OH/RscriptTesting/DDScatter/
ddplot=$HOME/devel/5OH/src/R/diffuse_discrete_scatter.R

# Loop samples and counts
n=${#SAMPLES}
for ((i=0;i<$n;i++)); do
    for count in "${MINCOUNTS[@]}"; do
        Rscript $ddplot "${SAMPLES[$i]}.intersect.tab" ${SAMPLES[$i]} ${DESCR[$i]} $count $webpath
    done
done