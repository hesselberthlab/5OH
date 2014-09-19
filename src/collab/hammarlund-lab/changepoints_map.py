#! /usr/bin/env python

''' changpoints_map: calc changepoints from genomecov map output. print '''

import sys
import ipdb
import operator

from numpy import mean, median, asarray, isnan
from segtools import ProgressBar

from rpy2 import robjects 
from rpy2.robjects.packages import importr

__version__ = '0.1'

# Import changepoint package from R
try:
    changepoint = importr('changepoint')
except ImportError:
    print >>sys.stderr, ">> need R: install.packages('changepoint')"
    sys.exit(1)

# minimum number of signals for changepoint calculation
CPT_MINSIGNALS = 4

# count cutoff requirement for input in changepoint, to eliminate intron noise
COUNT_CUTOFF = 10

def changepoints(pos_signal_mapping, neg_signal_mapping): #, verbose):
    parse_file(pos_signal_mapping)
    parse_file(neg_signal_mapping)


def parse_file(strand_signal_file):
    """Read in map file, parse collapsed mappings, pass to changepoints functions"""
    
    for line in open(strand_signal_file):
        fields = line.strip().split("\t")
        bed6fields = fields[:6]
        strand = fields[5]
        
        # Parse string of comma-separated count values to int list
        count_list = [int(x) for x in fields[6].split(",")]
        
        # Process gene if there are at least CPT_MINSIGNALS above COUNT_CUTOFF
        # TODO: put in minimum count cutoff or something?
        
        cl_array = asarray(count_list) # count_list as numpy array
        cl_above_cutoff = cl_array[cl_array >= COUNT_CUTOFF]
        
        if len(cl_above_cutoff) >= CPT_MINSIGNALS:
            process_changepoints(bed6fields, count_list, strand)
        else:
            print_zero_bed_line(bed6fields)


def process_changepoints(bed6fields, count_list, strand):
    '''calculate changepoints, calculate interval score, and print BED line'''

    interval_cpt = calc_changepoint(count_list)
    score = calc_interval_score(count_list,interval_cpt,strand)
    if interval_cpt and not isnan(score):
        print_cpt_bed_line(bed6fields,score,interval_cpt,strand)
    else:
        print_zero_bed_line(bed6fields)
    

def print_zero_bed_line(bed6fields):
    '''Print BED line with chrom, start, stop, gene name, zero score (no changepoint), strand'''

    chrom, start, stop, gene_name, count, strand = bed6fields
    fields = (chrom, start, stop, gene_name, 0, strand)
    print '\t'.join(map(str, fields))


def print_cpt_bed_line(bed6fields,score,cpt,strand):
    '''Print BED line with chrom, start, stop, gene name, cpt score, strand'''

    chrom, start, stop, gene_name, count, strand = bed6fields

    # Change stard and end coords to correspond to changepoint site
    if strand == "+":
        start = int(start) + cpt
        stop = start + 1
    elif strand == "-":
        stop = int(stop) - cpt
        start = stop - 1

    fields = (chrom, start, stop, gene_name, score, strand)
    print '\t'.join(map(str, fields))


def calc_interval_score(count_list, interval_cpt, strand):
    ''' calculates score for the interval based on a changepoint. score is
    the mean signal 3' of the cpt divided by the mean score 5' of the
    changepoint.
    
    Returns:
        score (float)
        '''
        
    # Splice count_list into signal left and right of changepoint
    left_cpt_cl = asarray(count_list[:interval_cpt])
    right_cpt_cl = asarray(count_list[interval_cpt:])
    
    # Ignore low signals (Consider better sort of penalty, etc)
    left_cpt_signal = mean(left_cpt_cl[left_cpt_cl >= COUNT_CUTOFF])
    right_cpt_signal = mean(right_cpt_cl[right_cpt_cl >= COUNT_CUTOFF])

    # Score takes strand into account, i.e. score = downstream / upstream 
    if strand == '+':
        score = right_cpt_signal / left_cpt_signal
    elif strand == '-':
        score = left_cpt_signal / right_cpt_signal

    return score

def calc_changepoint(count_list):
    '''Return first changepoint given a list of counts.
    
    Returns:
        cpoint (int)
        '''
    # Note: Currently only considers single (first) changepoint
    # Can change for multiple (or maybe last?) changepoint

    signals = robjects.IntVector(count_list)

    cpt_data = changepoint.cpt_mean(signals)
    cpoints = changepoint.cpts(cpt_data)

    # Get first calculated changepoint.  Cpoints is an Int Vector of changepoints.
    if not cpoints: return None
    cpoint = int(cpoints[0])

    return cpoint

def main():

    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(description=__doc__,
                            version=__version__,
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--pos-signal-mapping', help='pos mapping')
    parser.add_argument('-n', '--neg-signal-mapping', help='neg mapping')

    args = parser.parse_args()

    return changepoints(args.pos_signal_mapping,
                        args.neg_signal_mapping)

if __name__ == '__main__':
    sys.exit(main())

