#!/usr/bin/env bash

#BSUB -J master
#BSUB -e master.%J.err
#BSUB -o master.%J.out
#BSUB -q normal

<<DOC
master analysis loop for 5OH methods-paper pipeline
DOC

set -o nounset -o pipefail -o errexit -x

# Starting with a single sacCer1 assembly for troubleshooting;
# retaining option to expand to more assemblies later.
export ASSEMBLIES=("sacCer1")
#export ASSEMBLIES=("sacCer1" "sacCer2" "sacCer3")
export PIPELINE=$HOME/devel/5OH/src/pipeline
export CONFIG=$PIPELINE/config.sh

for assembly in ${ASSEMBLIES[@]}; do

    source $CONFIG

    # reassign assembly-specific variables
    export ASSEMBLY=$assembly
    # XXX DEBUG provides a directory extension ("common-debug"), or
    # nothing ("common")
    export RESULT=$HOME/projects/5OH/results/common$DEBUG/$assembly

    export BOWTIEIDX=/vol3/home/jhessel/ref/genomes/$assembly/$assembly
    export CHROM_SIZES=$HOME/ref/genomes/$assembly/$assembly.chrom.sizes
    export GTF=$HOME/ref/genomes/$assembly/sgdGene.$assembly.gtf
    export FASTA=$HOME/ref/genomes/$assembly/$assembly.fa
    export FULLGFF=$HOME/ref/genomes/sacCer1/sacCer1.fuller.bed

    # can change to dynamic creation of these files here or in 4_makeWindows for size control
    export MRNAWINDOWS="/vol2/home/speach/ref/genomes/sacCer1/sacCer1.mrna.20windows.bed"
    export UTRWINDOWS="/vol2/home/speach/ref/genomes/sacCer1/sacCer1.UTRs.2window.bed"
  
    job_array="[1-$NUM_SAMPLES]"

    # job names look like: align_sacCer1[1-8]
    bsub -J "align_$ASSEMBLY$job_array" \
        < $PIPELINE/1_align.sh

    # jobs are dependent among individual job indices
    bsub -J "coverage_$ASSEMBLY$job_array" \
        -w "done('align_$ASSEMBLY[*]')" \
        < $PIPELINE/2_coverage.sh 

    bsub -J "createRtab_$ASSEMBLY$job_array" \
        -w "done('coverage_$ASSEMBLY[*]')" \
        < $PIPELINE/3_createRtab.sh

    bsub -J "makeWindows_$ASSEMBLY$job_array" \
        -w "done('createRtab_$ASSEMBLY[*]')" \
        < $PIPELINE/4_makeWindows.sh

    # There will be some plots eventually.
    #bsub -J "plots_$ASSEMBLY$job_array" \
    #    -w "done('origin_anal_$ASSEMBLY[*]') && \
    #        done('nuc_freqs_$ASSEMBLY[*]') && \
    #        done('exp_anal_$ASSEMBLY[*]')" \
    #    < $PIPELINE/5_plots.sh

    bsub -J "tracklines_$ASSEMBLY" \
        < $PIPELINE/tracklines.sh

done

