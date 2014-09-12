#! /usr/bin/env python

''' changpoints: calc changepoints '''

import sys
import ipdb
import operator

from numpy import mean, median
from segtools import ProgressBar
from pybedtools import BedTool, cleanup

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

def changepoints(gene_bed, pos_signal_bedgraph, neg_signal_bedgraph, verbose):

    genes_bedtool = BedTool(gene_bed)
    pos_signal_bedtool = BedTool(pos_signal_bedgraph)
    neg_signal_bedtool = BedTool(neg_signal_bedgraph)

    if verbose:
        # make a progress bar
        num_genes = len(list(genes_bedtool))
        progress = ProgressBar(num_genes, label=">> calculating cpts: ")

    for interval in genes_bedtool:      

        if verbose: progress.next()

        intersect_data = signal_from_interval(interval, pos_signal_bedtool,
                                              neg_signal_bedtool, verbose)
        if not intersect_data: continue

        # get counts from the intersect data and convert to IntVector
        counts = ([int(score) for chrom, start, end, score in intersect_data])
        signals = robjects.IntVector(counts)

        interval_cpt = calc_changepoint(signals)

        if not interval_cpt: continue

        interval_score = calc_interval_score(counts, interval_cpt,
                                    interval.strand, verbose)

        # the cpt is an index in the coord list, look up the orginal
        # coords for the report
        coords = [(chrom, start, end) for chrom, start, end, score
                  in intersect_data]

        chrom, start, stop = coords[interval_cpt]

        # report BED6 format
        fields = (chrom, start, stop, interval.name, interval_score, interval.strand)
        print '\t'.join(map(str, fields))

    if verbose: progress.end()

def calc_interval_score(signal, interval_cpt, strand, verbose):
    ''' calculates score for the interval based on a changepoint. score is
    the mean signal 3' of the cpt divided by the mean score 5' of the
    changepoint.
    
    Returns:
        score (float)
        '''

    left_cpt_signal = mean(signal[:interval_cpt])
    right_cpt_signal = mean(signal[interval_cpt:])

    # score takes strand into account, i.e. score = downstream / upstream 
    if strand == '+':
        score = right_cpt_signal / left_cpt_signal
    elif strand == '-':
        score = left_cpt_signal / right_cpt_signal

    return score

def signal_from_interval(interval, pos_signal_bedtool, neg_signal_bedtool, verbose):
    ''' intersects signal from a given interval. converts signals to
    IntVector for changepoint calculation'''

    if interval.strand == '+':
        signal_bedtool = pos_signal_bedtool
    elif interval.strand == '-':
        signal_bedtool = neg_signal_bedtool

    interval_bedtool = BedTool([interval.fields[:3]])

    intersect = signal_bedtool.intersect(interval_bedtool, sorted=True)

    intersect_data = [i.fields for i in intersect]

    if len(intersect_data) < CPT_MINSIGNALS: return None

    # explicity cleanup files with pybedtools.cleanup()
    cleanup(remove_all=True)

    return intersect_data

def calc_changepoint(signals):
    """Return first changepoint given an IntVector of counts."""
    # Note: Currently only considers single (first) changepoint
    # Can change for multiple (or maybe last?) changepoint

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

    parser.add_argument('gene_bed', help='BED file of genes')
    parser.add_argument('-p', '--pos-signal-bedgraph', help='pos bedgraph signal')
    parser.add_argument('-n', '--neg-signal-bedgraph', help='neg bedgraph signal')
    parser.add_argument('--verbose', action='store_true',
                        help='be verbose',default=False)

    args = parser.parse_args()

    return changepoints(args.gene_bed, args.pos_signal_bedgraph,
                        args.neg_signal_bedgraph, args.verbose)

if __name__ == '__main__':
    sys.exit(main())

