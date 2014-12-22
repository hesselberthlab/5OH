#!/usr/bin/env bash

#BSUB -J NBD
#BSUB -e NBD.%J.%I.err
#BSUB -o NBD.%J.%I.out
#BSUB -q normal
#BSUB -n 6

<<DOC
Run Harigaya's NBD scripts on my daterz.
DOC

set -o nounset -o pipefail -o errexit -x

mappingdir=~/projects/5OH/results/methodspaper/proportions/
nbddir=~/projects/5OH/results/methodspaper/nbd/
sample=SP8.assembly-sacCer1.align-uniq
genemap=${mappingdir}${sample}.genemap.tab
nbdprefix=${nbddir}${sample}

python make_harigaya_input_file.py $genemap > ${nbdprefix}.genemap.hp
Rscript nbinom_fit.R ${nbdprefix}.genemap.hp ${nbdprefix}_nb.txt

perl adjust_p_value.pl --in=${nbdprefix}_nb.txt --out=${nbdprefix}_fdr0.01.txt --fdr=0.01
perl adjust_p_value.pl --in=${nbdprefix}_nb.txt --out=${nbdprefix}_fdr0.05.txt --fdr=0.05
