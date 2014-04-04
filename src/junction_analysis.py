#! /usr/bin/env python

# Assess base bias at junction
# Site of linkage is defined as "0"
# Negative index values are 5' (linker)
# Positive index values are 3' (5OH sequence)
# Intended input is umi-trimmed fastq file

import sys
from pysam import Samfile
import pandas as pd
import pandas.rpy.common as com
from collections import defaultdict

umi_sam=sys.argv[1]

umi_len = 8
index_depth = 8		# desired index depth for pairwise analysis
					# value must be <= umi_len
index_counts = {}
pair_counts = {}

def junction_analysis(umi_sam):
	"""Run all necessary functions for data analysis.

		Params:
			umi_sam (str): bam file after umitools trim & rmdup
	"""
	initiate_ds()
	for umi, seq in parse_samfile(umi_sam):
		compare_bases(umi, seq)

	write_r_dataframe(index_counts,"index_counts.txt")
	write_r_dataframe(pair_counts,"pair_counts.txt")
	write_r_dataframe(compress_pairs(pair_counts),"pair_counts_comp.txt")


def parse_samfile(samfilename):
	"""Function to parse UMI and read sequence from sam (bam) file
       
        Params:
            samfilename (str): sam filename
        Returns:
	        UMI (str), sequence (str)
	    
	    ### Currently expects umitools (trim, rmdup) corrected bam file.
	"""
	with Samfile(samfilename, 'rb') as sfile:
	    for read in sfile.fetch():
			umi = read.qname.split("_")[1]
			yield umi, read.seq

def initiate_ds():
	"""Create data structures"""
	for i in range(1,index_depth+1):
		index_counts[i]=defaultdict(int)
		index_counts[-i]=defaultdict(int)
		pair_counts["pair-"+str(i)]=defaultdict(int)

def compare_bases(umi, seq):
	"""Compare base frequencies individually and pairwise
	by index number
	"""
	for i in range(1,index_depth+1):
		# 5' end of linked sequence
		sequence_base = seq[umi_len+i-1]
		index_counts[i][sequence_base] += 1
			
		# 3' end of linker
		linker_base = umi[umi_len-i]
		index_counts[-i][linker_base] += 1

		# Bases equidistance from linkage
		pair_bases = linker_base + sequence_base
		pair_counts["pair-"+str(i)][pair_bases] += 1

def compress_pairs(pd):
	"""Combine pair counts with counts of reverse pairs, ie AT==TA
		pd ~ pair dictionary
	"""
	new_pd = {}

	for pair in pd:
		new_pd[pair] = defaultdict(int)

		for key in pd[pair]:
			keys = [key, key[::-1]]		# key and reverse key
			keys = sorted(keys)

			# Prevents doubling of palindromic vals 
			if keys[0]==keys[1]:
				new_pd[pair][keys[0]] += pd[pair][key]
		
			if keys[0] not in new_pd[pair]:
				new_pd[pair][keys[0]] += pd[pair][keys[0]] + \
										 pd[pair][keys[1]]

	return new_pd

def write_r_dataframe(count_dict,outputfile):
	"""Write R-compatible dataframe to output file"""
	df = pd.DataFrame.from_dict(count_dict)
	r_dataframe = com.convert_to_r_dataframe(df)

	f1 = open(outputfile,'w+')
	print >>f1, r_dataframe


junction_analysis(umi_sam)


