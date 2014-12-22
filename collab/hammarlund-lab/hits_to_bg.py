#! /usr/bin/env python

''' hits_to_bg: remove duplicate entries in converted BED file'''

# TODO: Make pretty

import sys

filename = sys.argv[1]
SCORE_CHANGE = int(sys.argv[2])

prev_start = 0

for line in open(filename):
    fields = line.strip().split('\t')
    start1 = fields[1]
    
    if start1 == prev_start1: continue # isoforms with same changepoint, incompatible w/ BW
    
    else:
        prev_start1 = start1        
        compare_conditions(fields)
        
        
def compare_conditions(fields):
    chrom, start1, stop1, gene, score1, strand = fields[:4]
    start2, stop2, gene, score2 = fields[6:]
    
    # type conversions
    score1 = float(score1)
    score2 = float(score2)
    
    score_diff = score2 - score1
    
    # same change point & test condition has greater cpt score by SCORE_CHANGE threshold
    # print changepoint and the score_diff
    if start1 == start2 and score_diff > SCORE_CHANGE:
        print "\t".join([chrom, start2, stop2, gene, str(score_diff), strand])

    # different change point, but test condition above SCORE_CHANGE threshold
    # print the changepoint, but the score_diff is not relevant, score only
    if start1 != start2 and score2 > SCORE_CHANGE:
        print "\t".join([chrom, start2, stop2, gene, str(score2), strand])
    
    # if neither criteria is met, gene is skipped
  
