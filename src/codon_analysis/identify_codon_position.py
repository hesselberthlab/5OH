#! /usr/bin/env python

'''identify_codon_position.py: 
    input: a bedtools intersect output file with additional columns
    output: bedfile where 'score' is the codon position of 5OH site'''

import sys

filename = sys.argv[1]

for line in open(filename):
    line = line.strip().split("\t")
    chrom, start, stop, gene, count = line[:5]
    strand = line[10]

    if strand == "+":
        gene_start = int(line[6])
        codon_position = (int(start) - gene_start)%3 + 1

    if strand == "-":
        gene_start = int(line[7])
        codon_position = (gene_start - int(stop))%3 + 1

    print "\t".join([chrom, start, stop, gene, str(codon_position),strand])
