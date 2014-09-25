#! /usr/bin/env python

''' dual_changepoints_map:
    - calc changepoint from genomecov map output
    - compare scores between two samples at given changepoint
    - output in BED format if over threshold'''

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
CPT_MINSIGNALS = 20

# count cutoff requirement for input in changepoint, to eliminate intron noise
COUNT_CUTOFF = 10

def dual_changepoints_map(signalfile_a, signalfile_b):
    """Read in map file, parse collapsed mappings, pass to changepoints functions"""
    
    bfile = open(signalfile_b)
    for line_a in open(signalfile_a):
        line_b = bfile.readline()
        
        fields_a = line_a.strip().split("\t")
        bed6fields = fields_a[:6]
        strand = fields_a[5]
        
        # Parse string of comma-separated count values to int list
        count_list_a = [int(x) for x in fields_a[6].split(",")]
        
        # Process gene if there are at least CPT_MINSIGNALS above COUNT_CUTOFF
        cl_array = asarray(count_list_a) # count_list as numpy array
        cl_above_cutoff = cl_array[cl_array >= COUNT_CUTOFF]
        
        if len(cl_above_cutoff) >= CPT_MINSIGNALS:
            fields_b = line_b.strip().split("\t")
            count_list_b = [int(x) for x in fields_b[6].split(",")]
            process_changepoints(bed6fields, count_list_a, count_list_b, strand)


def process_changepoints(bed6fields, count_list_a, count_list_b, strand):
    '''calculate changepoints, calculate interval score, and print BED line'''

    interval_cpt = calc_changepoint(count_list_a)
    score_a = calc_interval_score(count_list_a,interval_cpt,strand)
    if interval_cpt and not isnan(score_a):
        score_b = calc_interval_score(count_list_b,interval_cpt,strand)
        score_diff = score_a - score_b

        if score_diff >=2 and not isnan(score_diff):
            print_cpt_bed_line(bed6fields,score_diff,interval_cpt,strand)
    

def print_zero_bed_line(bed6fields):
    '''Print BED line with chrom, start, stop, gene name, zero score (no changepoint), strand'''

    chrom, start, stop, gene_name, count, strand = bed6fields
    fields = (chrom, start, stop, gene_name, 0, strand)
    print '\t'.join(map(str, fields))


def print_cpt_bed_line(bed6fields,score_diff,cpt,strand):
    '''Print BED line with chrom, start, stop, gene name, cpt score difference, strand'''

    chrom, start, stop, gene_name, count, strand = bed6fields

    # Change stard and end coords to correspond to changepoint site
    if strand == "+":
        start = int(start) + cpt
        stop = start + 1
    elif strand == "-":
        stop = int(stop) - cpt
        start = stop - 1

    fields = (chrom, start, stop, gene_name, score_diff, strand)
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
    left_cpt_signal = mean(left_cpt_cl[left_cpt_cl >= COUNT_CUTOFF][-60:])
    right_cpt_signal = mean(right_cpt_cl[right_cpt_cl >= COUNT_CUTOFF][:60])

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

    parser.add_argument('-a', '--signalfile-a', help='condition_a mapping file')
    parser.add_argument('-b', '--signalfile-b', help='condition_b mapping file')
    #parser.add_argument('-c', '--cpt-min-signal', help='changepoint minimum signal')
    #parser.add_argument('-x', '--count-cutoff', help='count_cutoff')

    args = parser.parse_args()

    return dual_changepoints_map(args.signalfile_a,
                        args.signalfile_b)
                        #args.cpt_min_signal,
                        #args.count_cutoff)

if __name__ == '__main__':
    sys.exit(main())

