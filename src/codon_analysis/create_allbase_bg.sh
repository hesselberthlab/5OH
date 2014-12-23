#! /usr/bin/env bash

# Create an allbase bed file with stranded data

sample=$1
bam=$sample.bam
allbaseposbg=$sample.allbase.pos.bg
allbasenegbg=$sample.allbase.neg.bg
CPMbasepos=$sample.allbase.CPM.pos.bg
CPMbaseneg=$sample.allbase.CPM.neg.bg

bedtools genomecov -ibam $bam -d -strand + \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $2 + 1, $3}' \
    > $allbaseposbg

bedtools genomecov -ibam $bam -d -strand - \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $2 + 1, $3}' \
    > $allbasenegbg

posreads=$(awk '{sum+=$4} END{print sum}' $allbaseposbg)
negreads=$(awk '{sum+=$4} END{print sum}' $allbasenegbg)

reads=$(($posreads + $negreads))

awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $allbaseposbg \
    > $CPMSposbg

awk -v total=$reads '{OFS="\t"} {printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4/total*1000000}' $allbasenegbg \
    > $CPMbaseng
