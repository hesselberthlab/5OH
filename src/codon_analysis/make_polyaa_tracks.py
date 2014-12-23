#! /usr/bin/env python

"""
Take the polyaa bed-like file and make some sensible motif tracks for the
UCSC browser
"""

import sys
import re

polyaabed = "polyall-RKED_genes.bed"

basic = ["R", "K", "[RK]", "[RKH]"]
acidic = ["D", "E", "[DE]"]
charged = ["[RD]", "[RE]", "[KE]", "[KD]", "[RKE]", "[RKD]", "[RKED]"]

motifdir = "/vol2/home/speach/ref/genomes/sacCer1/aa_motifs/"

basicbedfile = motifdir + "aa_basic.bed"
acidicbedfile = motifdir + "aa_acidic.bed"
chargedbedfile = motifdir + "aa_charged.bed"
otherbedfile = motifdir + "aa_etc.bed"

basicbed = open(basicbedfile, 'w+')
acidicbed = open(acidicbedfile, 'w+')
chargedbed = open(chargedbedfile, 'w+')
otherbed = open(otherbedfile, 'w+')

for line in open(polyaabed):
    chrom, start, stop, gene, motif, strand = line.strip().split("\t")
    motif_aminos = re.split('[0-9]', motif)[0] # parse aminos in motif 
    if motif_aminos.startswith("["):
        motif_search = "\[" + motif_aminos[1:-1] + "\]"
    else: motif_search = motif_aminos

    motif_length = int(re.split(motif_search, motif)[-1])

    if strand == "+":
        start = int(start) + 75
        stop = start + 3 * motif_length

    if strand == "-":
        stop = int(stop) - 75
        start = stop - 3 * motif_length

    bedline = "\t".join([chrom, str(start), str(stop), motif,
                        str(motif_length), strand]) + "\n"

    if motif_aminos in basic:
        basicbed.write(bedline)

    elif motif_aminos in acidic:
        acidicbed.write(bedline)

    elif motif_aminos in charged:
        chargedbed.write(bedline)

    else:
        otherbed.write(bedline)

basicbed.close()
acidicbed.close()
chargedbed.close()
otherbed.close()
