#! /usr/bin/env python

import sys
from collections import defaultdict
import pandas as pd
import numpy as np
import re


filename = sys.argv[1]
csv_filename = sys.argv[2]
toBin = False
if len(sys.argv) == 5:
    toBin = True
    NUM_BINS = sys.argv[4]

def create_matrices_dict(filename):
    '''Returns a dictionary of matrices keyed by poly-aa motif.  The matrix is a list of
        count vectors at each index, by gene.'''
    matrices_dict = defaultdict(list)
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
        countarray = map(float, countarray)
        
        # Reverse counts if - strand data
        if strand == "-":
            countarray.reverse()
        
        # Remove brackets from regex polyaa keys
        if polyaa_name.startswith("["):
            polyaa_name = polyaa_name.strip("[")
            polyaa_name = re.sub(']', '', polyaa_name)
        
        if toBin:
            bin_countarray = get_binned_countarray(countarray, NUM_BINS)
            matrices_dict[polyaa_name].append(bin_countarray)
        else:
            matrices_dict[polyaa_name].append(countarray)
        prev_countarray = countarray

    return matrices_dict

def get_binned_countarray(countarray, num_bins):
    bin_len = int(len(countarray)/num_bins)
    
    bin_countarray = []
    for i in range(num_bins):
        bin_mean = np.mean(countarray[(bin_len*i) : (bin_len*(i+1))])
        bin_countarray.append(bin_mean)
    
    return bin_countarray

def avg_polyaa_matrices(mdict, csv_filename):
    '''average indexes of each polyaa-motif, write to R-compatible csv'''
    if toBin:
        avgs_dict = {'bin':range(1,NUM_BINS+1)}
    else:
        avgs_dict = {'index':range(0,151)} # should not be hard coded
    
    print "\t".join(["mismatch", "num.motifs"])

    for polyaa, matrix in mdict.iteritems():
        df = pd.DataFrame(matrix)
        index_avgs = df.mean().tolist()
        avgs_dict[polyaa] = index_avgs
        #print "\t".join([polyaa.strip().split(":")[-1], str(len(matrix))])
        print "\t".join([polyaa.strip(), str(len(matrix))])

    avgs_df = pd.DataFrame(avgs_dict)
    avgs_df.to_csv(index=False, path_or_buf=csv_filename)


# make main at some point
mdict = create_matrices_dict(filename)
avg_polyaa_matrices(mdict, csv_filename)
