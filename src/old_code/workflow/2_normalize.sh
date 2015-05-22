#! /usr/bin/env bash

#BSUB -J norm[1-8]
#BSUB -o err-out/%J.%I.out
#BSUB -e err-out/%J.%I.err

set -o nounset -o pipefail -o errexit -x

source config.sh

reads=$(awk '{sum+=$4} END{print sum}' $BGFILE)

# TPM (Li, Dewey BMC Bioinformatics; Pachter CSH lecture) -> CPM
# count/total * 1000000
#
awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $BGFILE \
  > $NORM_BGFILE

awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $POSBGFILE \
  > $NORM_POSBGFILE

awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $NEGBGFILE \
  > $NORM_NEGBGFILE

# RPKM (reads per kilbase per million mapped reads) -- Tiny numbers :(
# count / ( (genome/1000) ) * (total/1000000) )

#reads=$(awk '{sum+=$4} END{print sum}' $BGFILE)
#awk -v total=$reads -v genome=$GENOMESIZE '{OFS="\t"} {print $1, $2, $3, $4/(genome/1000*total/1000000)*100}' $BGFILE \
#  > $NORM_BGFILE

# Normalized to rRNA - meh.
# RRNA_NORM=$(grep -P "chr12\t461758" $BGFILE | awk '{print $4}')
#awk -v rrna=$RRNA_NORM '{OFS="\t"} {printf($1 "\t" $2 "\t" $3 "\t")} {printf("%0.f\n", ($4/rrna*10000))}' $BGFILE \
#	> $NORM_BGFILE

# For an mRNA normalization...
# MRNA_NORM=$(grep -P "chr11\t279775" $BGFILE | awk '{print $4}')


