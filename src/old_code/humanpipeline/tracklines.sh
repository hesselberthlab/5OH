#! /usr/bin/env bash

#BSUB -J tracklines
#BSUB -o tracklines.%J.out
#BSUB -e tracklines.%J.err
#BSUB -q short

<<DOC
Generate tracklines for UCSC
DOC

set -o nounset -o pipefail -o errexit -x

source config.sh

urlbase="http://amc-sandbox.ucdenver.edu"
strands=("both" "pos" "neg")
umi_types=("removed" "UMIs_not_removed")

tracklinefile=$RESULT/tracklines.txt
if [[ ! -d $RESULT ]]; then
    mkdir -p $RESULT
fi

# remove old file if exists b/c we're cat'ing
if [[ -f $tracklinefile ]]; then
    rm -f $tracklinefile
fi


for sample_idx in ${!SAMPLES[@]}; do

    # provided by config.sh
    sample=${SAMPLES[$sample_idx]}
    color=${COLORS[$sample_idx]}
    descrip=${DESCRIPS[$sample_idx]}

    webroot="$urlbase/~speach/projects/5OH/results/$DATEDIR$DEBUG/$ASSEMBLY/$sample"

    bigwigdir="$webroot/bedgraphs"
    bigbeddir="$webroot/peaks"

    for align_mode in ${ALIGN_MODES[@]}; do
        for strand in ${strands[@]}; do

            track_color="color=$color"

            # bigWig tracks
            bigwig="$bigwigdir/${sample}.align.${align_mode}.strand.${strand}.counts.bw"
            bigwig_name="name='$sample $strand $align_mode'"
            bigwig_descrip="description='$descrip sample=$sample strand=$strand \
                            align.mode=$align_mode'"
            bigwig_height="maxHeightPixels=50:35:35"

            bigwig_track="track type=bigWig"
            bigwig_url="bigDataUrl=$bigwig"
            bigwig_vis="visibility=full"
            bigwig_trackline="$bigwig_track $bigwig_url $bigwig_name \
                              $bigwig_descrip $track_color \
                              $bigwig_vis $bigwig_height"

            # print the line
            echo $bigwig_trackline >> $tracklinefile

        done
    done
done
