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

def changepoints(gene_bed, signal_bedgraph, verbose):

    genes_bedtool = BedTool(gene_bed)
    signal_bedtool = BedTool(signal_bedgraph)

    if verbose:
        # make a progress bar
        num_genes = len(list(genes_bedtool))
        progress = ProgressBar(num_genes, label=">> calculating cpts: ")

    for interval in genes_bedtool:      

        if verbose: progress.next()

        result = signal_from_interval(interval, signal_bedtool, verbose)

        if not result: continue
        coords, signal = result

        interval_cpt = calc_changepoint(signal)

        if not interval_cpt: continue

        score = calc_interval_score(signal, interval_cpt, interval.strand, verbose)

        # the cpt is an index in the coord list, look up the orginal
        # coords for the report
        chrom, start, stop = coords[interval_cpt]

        # report BED6 format
        fields = (chrom, start, stop, interval.name, score, interval.strand)
        print '\t'.join(map(str, fields))

    if verbose: progress.end()

def calc_interval_score(signal, interval_cpt, strand, verbose):
    ''' calculates score for the interval based on a changepoint. score is
    the mean signal 3' of the cpt divided by the mean score 5' of the
    changepoint.'''

    left_cpt_signal = mean(signal[interval_cpt:])
    right_cpt_signal = mean(signal[:interval_cpt])

    # score takes strand into account, i.e. score = downstream / upstream 
    if strand == '+':
        score = right_cpt_signal / left_cpt_signal
    elif strand == '-':
        score = left_cpt_signal / right_cpt_signal

    return score

def signal_from_interval(interval, signal_bedtool, verbose):
    ''' doc '''

    interval_bedtool = BedTool([interval.fields[:3]])
    intersect = signal_bedtool.intersect(interval_bedtool, sorted=True)

    # Create list from count intersection
    counts = [int(datum.fields[3]) for datum in intersect]

    # need 4 values to calc cpt
    if len(counts) < 4: return None

    # Convert to Int Vector
    counts = robjects.IntVector(counts)

    # maintain list of coords for later
    coords = [(i.chrom, i.start, i.end) for i in intersect]

    # explicity cleanup files with pybedtools.cleanup()
    cleanup(remove_all=True)

    return (coords, counts)

def calc_changepoint(signal):
    """Return first changepoint given an IntVector of counts."""
    # Note: Currently only considers single (first) changepoint
    # Can change for multiple (or maybe last?) changepoint

    cpt_data = changepoint.cpt_meanvar(signal)
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
    parser.add_argument('signal_bedgraph', help='bedgraph signal')
    parser.add_argument('--verbose', action='store_true',
                        help='be verbose',default=False)

    args = parser.parse_args()

    return changepoints(args.gene_bed, args.signal_bedgraph, args.verbose)

if __name__ == '__main__':
    sys.exit(main())

