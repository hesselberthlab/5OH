#! /usr/bin/env python

#import rpy2.objects as robjects
#from rpy2.robjects.packages import importr
#mass = import("mass")
from collections import defaultdict
import sys

#genemapfile = "/vol2/home/speach/projects/5OH/results/methodspaper/proportions/SP8.assembly-sacCer1.align-uniq.genemap.tab"
genemapfile = sys.argv[1]

def print_tabbed_counts(gene, countarray):
    countlist = [gene]
    for count in countarray:
        count = str(count)
        countlist.append(count)
    
    print "\t".join(countlist)
        
for line in open(genemapfile):
    # parse line
    chrom, gene_start, gene_stop, gene, score, strand, counts = line.strip().split("\t")
    
    # Split map counts
    countarray = counts.split(",")
    
    if countarray == ["."]: continue # no data maps here
    
    else:
        print gene + "\t" + "\t".join(countarray)
