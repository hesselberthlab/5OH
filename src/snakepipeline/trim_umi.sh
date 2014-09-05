#!/usr/bin/env bash

#BSUB -J trim
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -n 6

<<DOC
Trim the UMI from the FASTQ.
Input: DATA + {sample}.fq.gz
Params: UMI
Output: DATA + {sample}.umi.fq.gz
DOC

set -o nounset -o pipefail -o errexit -x

# Read variables from command line
unprocessed_fastq=$1
UMI=$2
fastq=$3

umitools trim $unprocessed_fastq $UMI \
    | gzip -c \
    > $fastq