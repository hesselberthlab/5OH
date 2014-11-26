#!/usr/bin/env bash

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