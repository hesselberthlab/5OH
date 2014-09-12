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

    genes = BedTool(gene_bed)
    signal = BedTool(signal_bedgraph)

    if verbose:
        # make a progress bar
        num_genes = len(list(genes))
        progress = ProgressBar(num_genes, label=">> calculating cpts: ")

    for interval in genes:      

        if verbose: progress.next()

        signal = signal_from_interval(interval, signal_bedtool, verbose)

        if not signal: continue

        interval_cpt = calc_changepoint(signal)

        # XXX change calc_changepoint() method to return None if no
        # changepoints, should be no numerical logic here i.e.:
        # if not interval_cpt: continue
        #
        if interval_cpt < 2: continue

        score = calc_interval_score(signal, interval_cpt, interval.strand, verbose)

        # the cpt is an index in the coord list, look up the orginal
        # coords for the report
        chrom, start, stop = signal[interval_cpt].fields[:3]

        # report BED6 format
        fields = (chrom, start, stop, interval.name, score, interval.strand)
        print '\t'.join(map(str, fields))

    if verbose: progress.end()

def calc_interval_score(signal, interval_cpt, strand, verbose):
    ''' calculates score for the interval based on a changepoint. score is
    the mean signal 3' of the cpt divided by the mean score 5' of the
    changepoint.'''

    # XXX: use numpy.mean
    left_cpt_signal = signal[interval_cpt:]
    right_cpt_signal = signal[:interval_cpt]

    ipdb.set_trace()

    # score takes strand into account, i.e. score = upstream / downstream
    if strand == '+':
        score = right_cpt_mean / left_cpt_mean
    elif strand == '-':
        score = left_cpt_mean / right_cpt_mean

    return score

def signal_from_interval(interval, signal_bedtool, verbose):
    ''' doc '''

    interval_bedtool = BedTool([interval.fields[:3]])
    # XXX strandedness argument?
    #
    # there is no strand in bedgraph data, ideally pass in both pos
    # and neg data and determine which one you need from the BED
    # strand field
    intersect = signal_bedtool.intersect(interval_bedtool, sorted=True)

    # Create list from count intersection
    counts = [int(datum.fields[3]) for datum in intersect]

    # need 4 values to calc cpt
    if len(counts) < 4: return None

    # Convert to Int Vector
    counts = robjects.IntVector(counts)

    # explicity cleanup files with pybedtools.cleanup()
    cleanup(remove_all=True)

    return counts

def calc_changepoint(signal):
    """Return first changepoint given an IntVector of counts."""
    # Note: Currently only considers single (first) changepoint
    # Can change for multiple (or maybe last?) changepoint

    cpt_data = changepoint.cpt_meanvar(signal)
    # XXX: cpoints is a list?
    cpoints = changepoint.cpts(cpt_data)

    if len(cpoints) != 0:
        cpoint = int(changepoint.cpts(cpt_data)[0])  # parse changepoint values
    else:
        cpoint = 0    # ie, no changepoint found

    return cpoint

def main():

    from argparse import ArgumentParser, RawDescriptionHelpFormatter

    parser = ArgumentParser(description=__doc__,
                            version=__version__,
                            formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument('gene_bed', help='BED file of genes')
    parser.add_argument('signal_bedgraph', help='bedgraph signal')
    parser.add_argument('--verbose', action='store_true',
                        help='be verbose [default: %default]')

    args = parser.parse_args()

    return changepoints(args.gene_bed, args.signal_bedgraph, args.verbose)

if __name__ == '__main__':
    sys.exit(main())

