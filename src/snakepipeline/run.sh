#!/usr/bin/env bash

#BSUB -J snakemake
#BSUB -e err-out/snakemake.%J.%I.err
#BSUB -o err-out/snakemake.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 12

set -o nounset -o pipefail -o errexit -x

snakemake -j