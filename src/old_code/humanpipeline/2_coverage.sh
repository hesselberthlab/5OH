#!/usr/bin/env bash

#BSUB -J coverage[1-4]
#BSUB -e coverage.%J.%I.err
#BSUB -o coverage.%J.%I.out
#BSUB -q short

<<DOC
Calculate 5' end coverage from bam files; normalize to CPMs (counts per million reads)
DOC

set -o nounset -o pipefail -o errexit -x

source config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

results=$RESULT/$sample
bgresults=$RESULT/$sample/bedgraphs
if [[ ! -d $bgresults ]]; then
    mkdir -p $bgresults
fi

umi_types=("removed" "UMIs_not_removed")

for umi_type in ${umi_types[@]}; do

    for align_mode in ${ALIGN_MODES[@]}; do

        if [[ $umi_type == 'removed' ]]; then
            bam=$results/alignment/$sample.align.$align_mode.bam
            countsbg=$bgresults/$sample.align.$align_mode.strand.both.counts.bg
            countsposbg=$bgresults/$sample.align.$align_mode.strand.pos.counts.bg
            countsnegbg=$bgresults/$sample.align.$align_mode.strand.neg.counts.bg

            cpmsbg=$bgresults/$sample.align.$align_mode.strand.both.CPMs.bg
            cpmsposbg=$bgresults/$sample.align.$align_mode.strand.pos.CPMs.bg
            cpmsnegbg=$bgresults/$sample.align.$align_mode.strand.neg.CPMs.bg
        else
            bam=$results/alignment/$sample.$umi_type.align.$align_mode.bam
            countsbg=$bgresults/$sample.$umi_type.align.$align_mode.strand.both.counts.bg
            countsposbg=$bgresults/$sample.$umi_type.align.$align_mode.strand.pos.counts.bg
            countsnegbg=$bgresults/$sample.$umi_type.align.$align_mode.strand.neg.counts.bg

            cpmsbg=$bgresults/$sample.$umi_type.align.$align_mode.strand.both.CPMs.bg
            cpmsposbg=$bgresults/$sample.$umi_type.align.$align_mode.strand.pos.CPMs.bg
            cpmsnegbg=$bgresults/$sample.$umi_type.align.$align_mode.strand.neg.CPMs.bg
        fi

        if [[ ! -f $countsbg ]]; then
            bedtools genomecov -5 -bg -g $CHROM_SIZES -ibam $bam \
                > $countsbg
            bedtools genomecov -5 -bg -g $CHROM_SIZES -ibam $bam \
                -strand "+" > $countsposbg
            bedtools genomecov -5 -bg -g $CHROM_SIZES -ibam $bam \
                -strand "-" > $countsnegbg

            # sum total number of reads in aligned file
            reads=$(awk '{sum+=$4} END{print sum}' $countsbg)

            # normalize stranded bedgraphs to counts per million reads
            awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $countsbg \
                > $cpmsbg

            awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $countsposbg \
                > $cpmsposbg

            awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $countsnegbg \
                > $cpmsnegbg
        fi
    done
done

# create bigwigs
for bgfile in $(ls $bgresults/*.bg); do
    bwfile="$bgresults/$(basename $bgfile .bg).bw"
    bedGraphToBigWig $bgfile $CHROM_SIZES $bwfile
done

