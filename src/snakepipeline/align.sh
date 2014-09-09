#!/usr/bin/env bash

<<DOC
Align the reads!

Input: DATA + {sample}.umi.fq.gz
Params: BOWTIEIDX, ALIGN_ARGS
Output: DATA + {sample}.UMIs_not_removed.align.{align_mode}.bam
DOC

set -o nounset -o pipefail -o errexit -x

# Read variables from command line
if [[ $# == 3 ]]; then        # When alignment arguments are empty, only three variables are passed to script
    fastq=$1
    BOWTIEIDX=$2
    umibam=$3
    align_arg=""
fi

if [[ $# == 4 ]]; then
    fastq=$1
    BOWTIEIDX=$2
    align_arg=$3
    umibam=$4
fi


# Align to genome with bowtie, index with samtools
zcat $fastq \
    | bowtie $align_arg --sam $BOWTIEIDX -p 6 - \
    | samtools view -ShuF4 - \
    | samtools sort -o - $umibam.temp -m 8G \
        > $umibam
samtools index $umibam