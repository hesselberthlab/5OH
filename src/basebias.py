# Assess base bias at junction
# Site of linkage is defined as "0"
# Negative index values are 5' (linker)
# Positive index values are 3' (5OH sequence)
# Intended input is umi-trimmed fastq file

import sys
import pandas as pd
import pandas.rpy.common as com
from collections import defaultdict

fastq=sys.argv[1]

umi_len = 8
index_depth = 5		# desired index depth for pairwise analysis
					# value must be <= umi_len

index_counts = {}
pair_counts = {}

def initiate_ds():
	"""Create data structures"""
	for i in range(1,index_depth+1):
		index_counts[i]=defaultdict(int)
		index_counts[-i]=defaultdict(int)
		pair_counts["pair-"+str(i)]=defaultdict(int)

def compare_bases(fastq):
	"""Compare base frequencies individually and pairwise
	by index number
	"""
	line_num = 0

	for line in open(fastq):
		line_num += 1

		if line_num % 4 == 2:
			line = line.strip()
			for i in range(1,index_depth+1):
				# 5' end of linked sequence
				sequence_base = line[umi_len+i-1]
				index_counts[i][sequence_base] += 1
			
				# 3' end of linker
				linker_base = line[umi_len-i]
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


initiate_ds()
compare_bases(fastq)
write_r_dataframe(index_counts,"index_counts.txt")
write_r_dataframe(pair_counts,"pair_counts.txt")
write_r_dataframe(compress_pairs(pair_counts),"pair_counts_comp.txt")

