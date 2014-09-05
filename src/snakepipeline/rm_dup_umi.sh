#!/usr/bin/env bash

#BSUB -J rm_dup_umi
#BSUB -e rm_dup_umi.%J.%I.err
#BSUB -o rm_dup_umi.%J.%I.out
#BSUB -q normal
#BSUB -n 6

<<DOC
Remove duplicate UMIs
Input: RESULT + "{sample}.UMIs_not_removed.align." + ALIGN_MODE + ".bam"
Params: RESULT + "{sample}.umi-report.align." + ALIGN_MODE + ".bed.gz"
Output: RESULT + "{sample}.align." + ALIGN_MODE + ".bam"
DOC

set -o nounset -o pipefail -o errexit -x

# Read variables from command line
umibam=$1
umibed=$2
bam=$3

# Remove duplicate UMIs with umitools
umitools rmdup $umibam $bam | gzip -c > $umibed
samtools index $bam