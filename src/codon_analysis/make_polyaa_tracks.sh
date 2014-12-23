#! /usr/bin/env bash

<<DOC
Taking the 150bp motif stretches and making some sensible motif files for
the UCSC browser.
DOC

codedir="/vol2/home/speach/devel/5OH/src/codon_analysis"

python ${codedir}/make_polyaa_tracks.py

CHROMSIZES="/vol2/home/speach//ref/genomes/sacCer1/sacCer1.chrom.sizes"
motifdir="/vol2/home/speach/ref/genomes/sacCer1/aa_motifs"

aa_sets=("aa_basic" "aa_acidic" "aa_charged" "aa_etc")

for aa_set in ${aa_sets[@]}; do
    bed=${motifdir}/${aa_set}.bed
    cp $bed ${bed}.tmp
    sort -k1,1 -k2,2n ${bed}.tmp > $bed
    rm -rf ${bed}.tmp

    bb=${motifdir}/${aa_set}.bb
    bedToBigBed $bed $CHROMSIZES $bb
    cp $bb ~/public_html/ref/genomes/sacCer1/aa_motifs/.
done
