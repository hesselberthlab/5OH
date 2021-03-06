#!/usr/bin/env bash

<<DOC
Align the reads!

Input: DATA + {sample}.umi.fq.gz
Params: BOWTIEIDX, ALIGN_ARGS
Output: DATA + {sample}.UMIs_not_removed.align.{align_mode}.bam
DOC

set -o nounset -o pipefail -o errexit -x

# Read variables from command line
if [[ $# == 4 ]] ; then        # When alignment arguments are empty, only three variables are passed to script
    fastq=$1 
    umibam=$2
    alignstats=$3
    BOWTIEIDX=$4
    align_arg=""
fi

if [[ $# > 4 ]] ; then
    fastq=$1
    umibam=$2
    alignstats=$3
    BOWTIEIDX=$4
    align_arg=${@:5:$#}
fi


# Align to genome with bowtie, index with samtools
zcat $fastq \
    | bowtie $align_arg --un ${fastq}.unaligned --sam $BOWTIEIDX -p 6 - \
        2> $alignstats \
    | samtools view -ShuF4 - \
    | samtools sort -o - $umibam.temp -m 8G \
        > $umibam
samtools index $umibam
