#! /usr/bin/env bash

<<DOC
Compare changepoint bed files between two samples

Useage: compare_changepoints.sh $bed1 $bed2 $comparison.name
DOC

# Two BED files to compare
#bed1=RTCB_004_TGACCA_L002.changepoints.bed
#comparisonname=RTCB_TM.vs.RTCB

bed1=CNTRL_002_CGATGT_L002.changepoints.bed
comparisonname=RTCB_TM.vs.WT_TM

bed2=RTCB_TM_012_CTTGTA_L002.changepoints.bed
scorechange=2
hitsfile=${comparisonname}.hits.tab
tempbg=${comparisonname}.temp.bg
bedgraph=${comparisonname}.bg
bigwig=${comparisonname}.bw
CHROMSIZES=/vol2/home/speach/ref/genomes/ce10/ce10.chrom.sizes


# Print diff, sort by difference
paste $bed1 $bed2 \
  | awk 'BEGIN{OFS="\t"; print "chrom\tstart\tstop\tgene\tstrand\tcond1.cptscore\tcond2.cptscore\tdiff"} \
      {print $7, $8, $9, $10, $12, $5, $11, $11-$5}' \
  | sort -k8,8gr \
  > $hitsfile

# BED format (score is difference) -> BG format for score.diff above foldchange -> BW
paste $bed1 $bed2 \
  | awk 'BEGIN{OFS="\t"} {print $7, $8, $9, $10, $11-$5, $12}' \
  | awk -v change=$scorechange 'BEGIN{OFS="\t"} $5 > change {print $1, $2, $3, $5}' \
  | sort -k1,1 -k2,2n \
  > $tempbg
  
python rmdup_bg.py $tempbg > $bedgraph
rm -rf $tempbg
  
bedGraphToBigWig $bedgraph $CHROMSIZES $bigwig
cp $bigwig $HOME/public_html/projects/collab/hammarlund-lab/results/.

# For quickly pasting results for analysis...
#sort -k4,4gr $bedgraph \
#    | awk '{print $1 ":" $2}' \
#    | less