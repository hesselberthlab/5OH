#!/usr/bin/env bash

#BSUB -J master
#BSUB -e master.%J.err
#BSUB -o master.%J.out
#BSUB -q normal

<<DOC
master analysis loop for human 5OH pipeline
DOC

set -o nounset -o pipefail -o errexit -x

# can add additional assemblies if necessary
export ASSEMBLIES=("hg19")
export PIPELINE=$HOME/devel/5OH/src/pipeline
export CONFIG=config.sh

for assembly in ${ASSEMBLIES[@]}; do

    source $CONFIG

    # reassign assembly-specific variables
    export ASSEMBLY=$assembly
    # XXX DEBUG provides a directory extension ("common-debug"), or
    # nothing ("common")
    export RESULT=$HOME/projects/5OH/results/$DATEDIR$DEBUG/$assembly

    export BOWTIEIDX=/vol3/home/jhessel/ref/genomes/$assembly/$assembly
    export CHROM_SIZES=/vol3/home/jhessel/ref/genomes/$assembly/$assembly.chrom.sizes
    export GTF=/vol3/home/jhessel/ref/genomes/$assembly/sgdGene.$assembly.gtf
    export FASTA=/vol3/home/jhessel/ref/genomes/$assembly/$assembly.fa
  
    job_array="[1-$NUM_SAMPLES]"

    # job names look like: align_sacCer1[1-8]
    bsub -J "align_$ASSEMBLY$job_array" \
        < $PIPELINE/1_align.sh

    # jobs are dependent among individual job indices
    bsub -J "coverage_$ASSEMBLY$job_array" \
        -w "done('align_$ASSEMBLY[*]')" \
        < $PIPELINE/2_coverage.sh 

    # none of this is human compatible yet.  in due time, in due time.
    #bsub -J "createRtab_$ASSEMBLY$job_array" \
    #    -w "done('coverage_$ASSEMBLY[*]')" \
    #    < $PIPELINE/3_createRtab.sh

    #bsub -J "makeWindows_$ASSEMBLY$job_array" \
    #    -w "done('createRtab_$ASSEMBLY[*]')" \
    #    < $PIPELINE/4_makeWindows.sh

    #bsub -J "plots_$ASSEMBLY$job_array" \
    #    -w "done('origin_anal_$ASSEMBLY[*]') && \
    #        done('nuc_freqs_$ASSEMBLY[*]') && \
    #        done('exp_anal_$ASSEMBLY[*]')" \
    #    < $PIPELINE/5_plots.sh

    bsub -J "tracklines_$ASSEMBLY" \
        < $PIPELINE/tracklines.sh

done

