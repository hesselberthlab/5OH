#! /usr/bin/env bash

<<DOC
Compare changepoint mapping files between two samples

Useage: compare_changepoints.sh $map_a $map_b $comparisonname
DOC

# Read Arguments; Bed1 is condition we are interested in.
map_a=$1
map_b=$2
compname=$3

comparisonbed=${compname}.60.bed
tempbg=${compname}.temp.bg
bedgraph=${compname}.bg
bigwig=${compname}.bw
CHROMSIZES=/vol2/home/speach/ref/genomes/ce10/ce10.chrom.sizes
PIPELINE="/vol2/home/speach/devel/5OH/src/collab/hammarlund-lab/"

# Python script to skip duplicates and generate BED of hits above scorechange threshold
python ${PIPELINE}dual_changepoints_map.py -a $map_a -b $map_b > $comparisonbed

# Select for hits above scorechange threshold
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' $comparisonbed \
    | sort -k1,1 -k2,2n \
    > $tempbg
  
python ${PIPELINE}rmdup_bg.py $tempbg > $bedgraph
  
bedGraphToBigWig $bedgraph $CHROMSIZES $bigwig
cp $bigwig $HOME/public_html/projects/collab/hammarlund-lab/results/.