#!/usr/bin/env bash

# Project Descriptors specifics
# Project: Rep 01, WT/xrn1 +/- Tm, 5OH libraries
DATE=20140107
PROJECTDATA="$HOME/projects/5OH/data/${DATE}/"
#PROJECTDIR="projects/5OH/results/${DATE}"
PROJECTDIR="projects/5OH/results/20140415_rep1"
BARCODES=(GCCTAA_S1 TGGTCA_S2 CACTGT_S3 ATTGGC_S4 TCAAGT_S5 CTGATC_S6 AAGCTA_S7 GTAGCC_S8)
DESCRS=("delXRN1 +Tm" "delXRN1 DMSO" "WT +Tm" "WT DMSO" "delXRN1 +Tm +SAP" "delXRN1 DMSO +SAP" "WT +Tm +SAP" "WT DMSO +SAP")

# Reference information
BOWTIEINDEX="/vol3/home/jhessel/ref/genomes/sacCer1/sacCer1"
CHROMSIZES="/vol3/home/jhessel/ref/genomes/sacCer1/sacCer1.chrom.sizes"
FASTA="/vol3/home/jhessel/ref/genomes/sacCer1/sacCer1.fa"
GENOMESIZE=12157105
MRNABED="/vol2/home/speach/ref/genomes/sacCer1/lee-segments.sorted.bed"
FULLGFF="/vol2/home/speach/ref/genomes/sacCer1/sacCer1.fuller.bed"

# Variables for UCSC Mapping
WEBFOLDER=$HOME/public_html/
SANDBOX="http://amc-sandbox.ucdenver.edu/~speach/"
RED="215,25,28"
ORANGE="253,174,97"
GREEN="0,100,0"
BLUE="43,131,186"
COLORS=($RED $ORANGE $GREEN $BLUE $RED $ORANGE $GREEN $BLUE)

# LSB indexing
numBarcodes=${#BARCODES[@]}
BARCODE=${BARCODES[$(($LSB_JOBINDEX-1))]}
COLOR=${COLORS[$(($LSB_JOBINDEX-1))]}
DESCR=${DESCRS[$(($LSB_JOBINDEX-1))]}

# BG, BW, TAB files
BGFILE=$BARCODE.bg
POSBGFILE=$BARCODE.pos.bg
NEGBGFILE=$BARCODE.neg.bg
NORM_BGFILE=$BARCODE.norm.bg
NORM_POSBGFILE=$BARCODE.norm.pos.bg
NORM_NEGBGFILE=$BARCODE.norm.neg.bg
INTERSECT=$BARCODE.intersect.tab


