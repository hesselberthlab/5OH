#!/usr/bin/env bash

<<DOC
Calculate 5' end coverage from bam files,
Normalize to CPMs (counts per million reads)

Input: RESULT + "{sample}.align." + ALIGN_MODE + ".bam"
Params: CHROM_SIZES, RESULT
Output: RESULT + "{sample}.align." + ALIGN_MODE + ".strand.all.CPMs.bg",
RESULT + "{sample}.align." + ALIGN_MODE + ".strand.pos.CPMs.bg",
RESULT + "{sample}.align." + ALIGN_MODE + ".strand.neg.CPMs.bg"
DOC

set -o nounset -o pipefail -o errexit -x

# Read variables from command line
bam=$1
CHROM_SIZES=$2
RESULT=$3
WEBFOLDER=$4
PROJECTDIR=$5
cov_args=""
if [[ $# > 5 ]] ; then
    cov_args=${@:6:$#}
fi
bamprefix=$(basename $bam .bam)
fileprefix=$RESULT/$bamprefix

# Create Bedgraph filenames
countsbg="$fileprefix.strand.all.counts.bg"
countsposbg="$fileprefix.strand.pos.counts.bg"
countsnegbg="$fileprefix.strand.neg.counts.bg"
cpmsbg="$fileprefix.strand.all.CPMs.bg"
cpmsposbg="$fileprefix.strand.pos.CPMs.bg"
cpmsnegbg="$fileprefix.strand.neg.CPMs.bg"


# Calculate coverage using bedtools genomecov
bedtools genomecov $cov_args $CHROM_SIZES -ibam $bam \
    > $countsbg
bedtools genomecov $cov_args $CHROM_SIZES -ibam $bam \
    -strand "+" > $countsposbg
bedtools genomecov $cov_args $CHROM_SIZES -ibam $bam \
    -strand "-" > $countsnegbg

# Sum total number of reads in aligned file using awk
reads=$(awk '{sum+=$4} END{print sum}' $countsbg)

# Normalize stranded bedgraphs to counts per million reads
awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $countsbg \
    > $cpmsbg

awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $countsposbg \
    > $cpmsposbg

awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $countsnegbg \
   > $cpmsnegbg

# Create bigwigs from all bedgraphs and copy to webfolder
if [[ ! -d ${WEBFOLDER}${PROJECTDIR} ]] ; then
    mkdir -p ${WEBFOLDER}${PROJECTDIR}/
fi

for bgfile in $(ls $fileprefix*.bg); do
    bwfile="$RESULT/$(basename $bgfile .bg).bw"
    bedGraphToBigWig $bgfile $CHROM_SIZES $bwfile
    cp $bwfile ${WEBFOLDER}${PROJECTDIR}/.
done
