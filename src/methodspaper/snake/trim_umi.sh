#!/usr/bin/env bash

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

# Trim UMI from fastq using umitoools trim
umitools trim $unprocessed_fastq $UMI \
    | gzip -c \
    > $fastq