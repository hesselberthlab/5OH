#! /usr/bin/env python

# Input: FASTA-like file in following format:
# name;chr:start-stop;[+|-] \t sequence
# (ie, ~speach/projects/5OH/results/20140304/sacCer1MRNA_s.fa)
#
# Output: Prints bed-formatted text with score field as highest probability
# of base pairing as computed by RNAplfold from ViennaRNA

import sys
from collections import defaultdict
from subprocess import Popen, PIPE

def generate_pair_probs(seq):
	"""Generates pair probabilities for given sequence
		Params:
			seq (str): RNA sequence
		Return:
			prob_dict (dict): 	key is base index
								value is highest prob of base pairing (local)
	"""
	p = Popen("RNAplfold", stdin=PIPE)		# part of ViennaRNA
	p.communicate(seq+'\n@')				# input sequence and exit prog

	table_start = False
	prob_dict = defaultdict(int)

	# process PostScript output file from RNAplfold and grab probabilities
	for line in open("plfold_dp.ps"):
		line = line.strip()
		if line == "showpage": break		# end of table
		if line == "drawgrid_turn":			# beginning of table
			table_start = True
			continue
		if table_start:
			i, j, prob = line.split(" ")[:3]
			if prob > prob_dict[i]:			# get greatest prob; (!!!not ideal soln)
				prob_dict[i] = prob

	return prob_dict

def make_bedfile(prob_table, label):
	"""Create a bedfile output from basepair probabilities and label
		Params:
			prob_table (dict) - dictionary of highest prob of base-pairing for given index
			label (str) - name;chr:start:stop;[+|-]
		Return:
			None
	"""
	name, coords, strand = label.split(";")
	chrom, startstop = coords.split(":")
	start, stop = startstop.split("-")
	start = int(start)
	
	for key in sorted(prob_table.keys()):
		i = int(key)
		print "\t".join([chrom, str(i+start-1), str(i+start), prob_table[key], name, strand])

def main(args=sys.argv[1:]):
	for line in open(args[0]):
		label, sequence = line.strip().split("\t")
		prob_dict = generate_pair_probs(sequence)
		make_bedfile(prob_dict, label)

if __name__ == '__main__':
	sys.exit(main()) 



