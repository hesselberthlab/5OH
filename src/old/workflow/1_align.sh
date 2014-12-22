#! /usr/bin/env bash

#BSUB -J align[1-8]
#BSUB -o err-out/%J.%I.out
#BSUB -e err-out/%J.%I.err

source config.sh

set -o nounset -o pipefail -o errexit -x

BAMPREFIX=$BARCODE
bam="$BAMPREFIX.bam"
umibam="$BAMPREFIX.umi.bam"
STDIN="-"

UMI="NNNNNNNN"
unproc_fastq="$PROJECTDATA/$BARCODE""_L001_R1_001.fastq.gz"
umi_fastq="$BARCODE.umi.fq.gz"

if [[ ! -f $umi_fastq ]]; then
    umitools trim $unproc_fastq $UMI | gzip -c > $umi_fastq
fi

if [[ ! -f $bam ]]; then
    zcat $umi_fastq \
        | bowtie -m 2 --sam $BOWTIEINDEX $STDIN \
        | samtools view -F 4 -bS $STDIN > "$BAMPREFIX.unsorted.bam"

    samtools sort "$BAMPREFIX.unsorted.bam" $BAMPREFIX
    samtools index $bam

    # remove unsorted bam file
    rm -f "$BAMPREFIX.unsorted.bam"
fi

if [[ ! -f $umibam ]]; then
    umitools rmdup $bam $umibam
    samtools index $umibam
fi

BEDGRAPH="$BAMPREFIX.bg"
BEDGRAPHPOS="$BAMPREFIX.pos.bg"
BEDGRAPHNEG="$BAMPREFIX.neg.bg"

bedtools genomecov -ibam $umibam \
    -5 -bg -g $CHROMSIZES \
    > $BEDGRAPH
bedtools genomecov -ibam $umibam \
    -5 -bg -g $CHROMSIZES -strand "+" \
    > $BEDGRAPHPOS
bedtools genomecov -ibam $umibam  \
    -5 -bg -g $CHROMSIZES -strand "-" \
    > $BEDGRAPHNEG
