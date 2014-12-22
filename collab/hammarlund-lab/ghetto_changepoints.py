#! /usr/bin/env python

''' ghetto_changpoints: calc changepoints, but not using pybedtools damnit. '''

import sys
import ipdb
import operator

from numpy import mean, median
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

def changepoints(pos_signal_intersect, neg_signal_intersect): #, verbose):
    parse_file(pos_signal_intersect)
    parse_file(neg_signal_intersect)


def parse_file(strand_signal_file):
    """Read in intersection file and create subbeds"""

    prev_gene = "none"
    subbeds = []
    count_list = []

    for line in open(strand_signal_file):
        fields = line.strip().split("\t")
        gene_name = fields[3]
        count = int(fields[4])
        strand = fields[5]

        if prev_gene == "none": # Initiation of the beautiful cycle
              prev_gene = gene_name
              subbeds.append(fields)
              count_list.append(count)

        elif gene_name == prev_gene: # Still on same gene
              subbeds.append(fields)
              count_list.append(count)

        elif gene_name != prev_gene: # Moved to a new gene
              if len(count_list) > CPT_MINSIGNALS:
                  process_changepoints(subbeds, count_list, strand)

              prev_gene = gene_name
              subbeds = []
              count_list = []
              subbeds.append(fields)
              count_list.append(count)

        # end condition for last line?  meh.

    #if verbose: progress.end()

def process_changepoints(subbeds, count_list, strand):
    '''calculate changepoints, calculate interval score, and print BED line'''

    interval_cpt = calc_changepoint(count_list)
    score = calc_interval_score(count_list,interval_cpt, strand)
    if interval_cpt: print_cpt_bed_line(subbeds,interval_cpt,score)

def print_cpt_bed_line(subbeds,interval_cpt,score):
    '''Print BED line with chrom, start, stop, gene name, cpt score, strand'''

    chrom, start, stop, gene_name, count, strand = subbeds[interval_cpt] # +/- 1?
    fields = (chrom, start, stop, gene_name, score, strand)
    print '\t'.join(map(str, fields))

def calc_interval_score(count_list, interval_cpt, strand):
    ''' calculates score for the interval based on a changepoint. score is
    the mean signal 3' of the cpt divided by the mean score 5' of the
    changepoint.
    
    Returns:
        score (float)
        '''
    left_cpt_signal = mean(count_list[:interval_cpt])
    right_cpt_signal = mean(count_list[interval_cpt:])

    # score takes strand into account, i.e. score = downstream / upstream 
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

    cpt_data = changepoint.cpt_meanvar(signals)
    cpoints = changepoint.cpts(cpt_data)

    if not cpoints: return None

    # first cpt in the list is the calculated changepoint. the last value
    # is the end of the signal
    cpoint = int(cpoints[0])

    return cpoint

def main():

    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(description=__doc__,
                            version=__version__,
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--pos-signal-intersect', help='pos mapping')
    parser.add_argument('-n', '--neg-signal-intersect', help='neg mapping')
    #parser.add_argument('--verbose', action='store_true',
    #                    help='be verbose',default=False)

    args = parser.parse_args()

    return changepoints(args.pos_signal_intersect,
                        args.neg_signal_intersect) #, args.verbose)

if __name__ == '__main__':
    sys.exit(main())

