#! /usr/bin/env python

''' hits_to_bg: remove duplicate entries in converted BED file'''

import sys
import ipdb

__version__ = '0.1'

def doublebed_to_hitsbed(filename, score_change):
    prev_start1 = 0
    
    for line in open(filename):
        fields = line.strip().split('\t')
        start1 = fields[1]
    
        if start1 == prev_start1: continue # isoforms with same changepoint, incompatible w/ BW
    
        else:
            prev_start1 = start1        
            compare_conditions(fields, score_change)
        
def compare_conditions(fields, score_change):
    chrom, start1, stop1, gene, score1, strand = fields[:6]
    start2, stop2, gene, score2 = fields[7:11]
    
    # type conversions
    score1 = float(score1)
    score2 = float(score2)
    
    score_diff = score2 - score1
    
    # same change point & test condition has greater cpt score by score_change threshold
    # print changepoint and the score_diff
    if start1 == start2 and score_diff > score_change:
        print "\t".join([chrom, start2, stop2, gene, str(score_diff), strand])

    # different change point, but test condition above score_change threshold
    if start1 != start2 and score2 > score_change:
        print "\t".join([chrom, start2, stop2, "***" + gene, str(score_diff), strand])
    
    # if neither criteria is met, gene is skipped
    
def main():

    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(description=__doc__,
                            version=__version__,
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--filename', help='doublebed filename')
    parser.add_argument('-s', '--score-change', help='score change threshold')

    args = parser.parse_args()

    return doublebed_to_hitsbed(args.filename,
                        int(args.score_change))

if __name__ == '__main__':
    sys.exit(main())
  
