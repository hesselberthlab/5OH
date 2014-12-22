#! /usr/bin/env python

''' rmdup_bg: remove duplicate entries in converted bedgraph file'''

import sys

filename = sys.argv[1]

prev_start = 0
firstline = True

for line in open(filename):
    if firstline:
        firstline = False
        continue
    chrom, start, stop, gene, score, strand = line.strip().split('\t')[:6]
    if start == prev_start: continue
    else:
        prev_start = start
        print "\t".join([chrom, start, stop, gene, score, strand])
  
