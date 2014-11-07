#! /usr/bin/env python

import sys
from collections import defaultdict
import pandas as pd
import re

filename = sys.argv[1]
csv_filename = sys.argv[2]

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
        countarray = map(int, countarray)
        
        # Reverse counts if - strand data
        if strand == "-":
            countarray.reverse()
        
        # Remove brackets from regex polyaa keys
        if polyaa_name.startswith("["):
            polyaa_name = polyaa_name.strip("[")
            polyaa_name = re.sub(']', '', polyaa_name)
            
        matrices_dict[polyaa_name].append(countarray)
        prev_countarray = countarray

    return matrices_dict


def avg_polyaa_matrices(mdict, csv_filename):
    '''average indexes of each polyaa-motif, write to R-compatible csv'''
    avgs_dict = {'index':range(0,101)}
    
    for polyaa, matrix in mdict.iteritems():
        df = pd.DataFrame(matrix)
        index_avgs = df.mean().tolist()
        avgs_dict[polyaa] = index_avgs
        r_df = pd.DataFrame({polyaa:index_avgs})
        
    avgs_df = pd.DataFrame(avgs_dict)
    avgs_df.to_csv(index=False, path_or_buf=csv_filename)


# make main at some point
mdict = create_matrices_dict(filename)
avg_polyaa_matrices(mdict, csv_filename)
