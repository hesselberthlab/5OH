#! /usr/bin/env bash

#BSUB -J feature_dens[1-4]
#BSUB -o feature_dens.%J.%I.out
#BSUB -e feature_dens.%J.%I.err

set -o nounset -o pipefail -o errexit -x

PROJECT=$HOME/devel/5OH
BIN=$PROJECT/src
RESULT=$HOME/projects/5OH/results/20140701
CHROMSIZE=$HOME/ref/genomes/hg19/hg19.chrom.sizes

ANNOT_TYPES=("polyAdb"
             'refGene.tss'
             'refGene.start-codons'
             'refGene.stop-codons')

ANNOT_FILES=("$PROJECT/data/annotations/hg19/polyaDb.size1.bed.gz"
             "$PROJECT/data/annotations/hg19/refGene.tss.bed.gz"
             "$PROJECT/data/annotations/hg19/refGene.start.codons.bed.gz"
             "$PROJECT/data/annotations/hg19/refGene.stop.codons.bed.gz")

samples=(scr sh1 sh2 sh3)
sample=${samples[$(($LSB_JOBINDEX - 1))]}

libtypes=(PolyA 5OH-Seq)

for libtype in ${libtypes[@]}; do

    posbedgraph="$RESULT/$libtype/$libtype""_$sample.pos.norm.bg"
    negbedgraph="$RESULT/$libtype/$libtype""_$sample.neg.norm.bg"

    for annot_idx in ${!ANNOT_TYPES[@]}; do
 
        annot_bed=${ANNOT_FILES[$annot_idx]}
        annot_type=${ANNOT_TYPES[$annot_idx]}

        result="$libtype.$sample.$annot_type.tab"

        args="--window-resolution=25 --feature-label=$annot_type
              --sample-label=$sample --library-type=$libtype"

        if [[ $libtype == "PolyA" ]]; then
                args="$args --invert-strand"
        fi

        cmd="python $BIN/feature_density.py \
            -p $posbedgraph \
            -n $negbedgraph \
            $annot_bed \
            $CHROMSIZE \
            $args \
            --verbose \
            > $result"
       
        jobname="feature_dens.$libtype.$sample.$annot_type"
        bsub -o %J.out -e %J.err -J $jobname $cmd

    done
done

