#! /usr/bin/env python

import sys

filename = sys.argv[1]

prev_counts = ""

for line in open(filename):
    line = line.strip().split("\t")
    polyaa_name = line[4]
    strand = line[5]
    counts = line[6]
    
    if counts == prev_counts:
        continue
    if counts == '.':
        continue
    
    countarray = counts.split(",")
    
    if strand == "-":
        countarray.reverse()

    tab_counts = "\t".join(countarray)
    
    if polyaa_name == "E5":
        print tab_counts
    
    prev_countarray = countarray
