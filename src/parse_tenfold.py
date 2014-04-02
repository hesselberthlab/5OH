#! /usr/bin/env python

# Generate BED-like file where "score" is tenfold output

import re
import sys

def parse_tenfold(tabfile):
	"""Convert tenfold.py output to BED output
	"score" is tenfold output:
	1 = double-stranded
	2 = single-standed
	"""
	for line in open(tabfile):
		label, secondary = line.strip().split("\t")

		split_label = re.findall(r"[\w']+", label)
		name = split_label[0]
		chrom, start, stop = split_label[-3:]

		start = int(start)
		stop = int(stop)

		for i in range(len(secondary)):
			print "\t".join([chrom, str(start+i), \
						str(start+i+1), name, secondary[i]])

parse_tenfold(sys.argv[1])
		
