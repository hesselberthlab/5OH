#! /usr/bin/env python

#import rpy2.objects as robjects
#from rpy2.robjects.packages import importr
#mass = import("mass")
from collections import defaultdict
import sys

#genemapfile = "/vol2/home/speach/projects/5OH/results/methodspaper/proportions/SP8.assembly-sacCer1.align-uniq.genemap.tab"
genemapfile = sys.argv[1]

def print_count_sums(countarray, gene): 
    countdict = defaultdict(int)
    for count in countarray:
        countdict[count] += 1

    for key, value in countdict.iteritems():
        print "\t".join([gene, str(key), str(value)])

def print_all_counts(chrom, gene_start, gene, countarray):
    gene_start = int(gene_start) - 1 # will increment for each site
    for count in countarray:
        gene_start += 1
        site = chrom + ":" + str(gene_start)
        print "\t".join([gene, site, str(count)])
        
        

#print "\t".join(["gene","count.level","num.counts"])
print "\t".join(["gene", "site", "count.level"])
for line in open(genemapfile):
    # parse line
    chrom, gene_start, gene_stop, gene, score, strand, counts = line.strip().split("\t")
    gene_start = int(gene_start)
    gene_stop = int(gene_stop)
    
    # Split map counts
    countarray = counts.split(",")
    
    if countarray == ["."]: continue # no data maps here
    
    else:
        countarray = map(int, countarray)
        print_all_counts(chrom, gene_start, gene, countarray)
