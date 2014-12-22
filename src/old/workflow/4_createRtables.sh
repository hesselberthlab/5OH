#! /usr/bin/env bash

#BSUB -J rtab[1-8]
#BSUB -o err-out/%J.%I.out
#BSUB -e err-out/%J.%I.err

set -o nounset -o pipefail -o errexit -x

source config.sh

bedtools intersect -a $POSBGFILE -b $FULLGFF -wao \
	| awk 'BEGIN{OFS="\t"; print "count\tcat\tgene"} $10=="+" {split($8, a, ":"); print $4, a[1], a[2]}' \
	> $INTERSECT

bedtools intersect -a $NEGBGFILE -b $FULLGFF -wao \
	| awk 'BEGIN{OFS="\t"} $10=="-" {split($8, a, ":"); print $4, a[1], a[2]}' \
	>> $INTERSECT
