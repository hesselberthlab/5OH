#! /usr/bin/env bash

#BSUB -J plots[1-8]
#BSUB -o plots.%J.%I.out
#BSUB -e plots.%J.%I.err

<<DOC
Run several Rscripts to generate plots of 5OH data
DOC

# WORK IN PROGRESS

set -o nounset -o pipefail -o errexit -x

source $CONFIG

# Parse sample-specific data from config file
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
samplename=${NAMES[$(($LSB_JOBINDEX - 1))]}
sampledescr=${DESCRIPS[$(($LSB_JOBINDEX - 1))]}

rtables=$RESULT/$sample/rtables
rplots=$RESULT/$sample/rplots
if [[ ! -d $rplots ]]; then
    mkdir -p $rplots
fi

samplefile=${rtables}/${sample}.

################################
# DIFFUSE/DISCRETE SCATTER PLOTS
#
# minimum number of counts
mincounts=(5 10 25 50 100)

for rtabfile in $(ls $rtables/*intersect.tab); do
    for count in "${mincounts[@]}"; do
      Rscript "${RSCRIPTS}/diffuse_discrete_scatter.R" $rtabfile $sample $sampledescr $count $rplots
done
###############################


###############################
# IDENTIFY SITE CHANGES
#
# which comparisons to make amongst samples?
# rep vs. rep; trmt vs. cntrl; etc
#
# [Rscript ...]
#
###############################
