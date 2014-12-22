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

resultsdir=$(dirname $umibam)
unalignedfq=${resultsdir}/$(basename $umibam .bam).unaligned.fq
trimmedunaligned=${resultsdir}/$(basename $umibam .bam).unaligned.trimmed.fq.gz
alignstats2=${resultsdir}/$(basename $alignstats .txt).2ndalignment.txt
umibam2=${resultsdir}/$(basename $umibam .bam).2ndalignment.bam

# Align to genome with bowtie, index with samtools
zcat $fastq \
    | bowtie $align_arg --un $unalignedfq --sam $BOWTIEIDX -p 6 - \
        2> $alignstats \
    | samtools view -ShuF4 - \
    | samtools sort -o - $umibam.temp -m 8G \
        > $umibam
samtools index $umibam

# Trim unaligned reads by 12 bases (42->30 and 3' end of read)
nbases="NNNNNNNNNNNN"
umitools trim $unalignedfq --end 3 $nbases \
    | gzip -c \
    > $trimmedunaligned

# Run 2nd alignment with the previously unaligned data
zcat $trimmedunaligned \
    | bowtie $align_arg --sam $BOWTIEIDX -p 6 - \
        2> $alignstats2 \
    | samtools view -ShuF4 - \
    | samtools sort -o - $umibam2.temp -m 8G \
        > $umibam2
samtools index $umibam2
