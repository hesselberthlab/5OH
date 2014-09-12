#! /usr/bin/env python

''' changpoints: calc changepoints '''

import sys
import pdb
import rpy2.robjects as robjects
import operator
import pybedtools

from segtools import ProgressBar
from pybedtools import BedTool
from rpy2.robjects.packages import importr

__version__ = '0.1'

# Import Changepoint package from R
try:
    changepoint = importr('changepoint')
except ImportError:
    print >>sys.stderr, "install R changepoint package"
    sys.exit(1)

def changepoints(gene_bed, signal_bedgraph, verbose):

    genes = BedTool(gene_bed)
    signal = BedTool(signal_bedgraph)

    changepoint_scores = {}

    if verbose:
        num_genes = len([region for region in genes])
        progress = ProgressBar(num_genes,
                              label=">> calculating cpts: ")

    for region in genes:      # loop through each gene
        gene_name = region.fields[3]
        region_bedtool = BedTool([region.fields[:3]])

        if verbose: progress.next()

        # XXX strandedness argument?
        signal_data = signal.intersect(region_bedtool, sorted=True)

        # Create list from count intersection
        count_list = [int(datum.fields[3]) for datum in signal_data]
        if len(count_list) < 4:               # need 4 values to calc changepoint
            pybedtools.cleanup(remove_all=True)
            continue

        # Convert to Int Vector
        count_vector = robjects.IntVector(count_list)

        region_cpt = calc_changepoint(count_vector)

        if region_cpt < 2:  # if no changepoint found
            pybedtools.cleanup(remove_all=True)
            continue

        # avg signal counts to right of changepoint
        rcpt_avg = sum(count_vector[region_cpt:])/float(len(count_vector[region_cpt:]))

        # avg signal counts to left of changepoint
        lcpt_avg = sum(count_vector[:region_cpt])/float(len(count_vector[:region_cpt]))

        # score is the ratio of avg. signal to right and left of changepoint
        score = rcpt_avg/lcpt_avg

        # Store in dict keyed by gene name
        # Consider other options, ie full region entry?
        # changepoint_scores[gene_name] = score

        region_cpt_coords = signal_data[region_cpt].fields[:3]

        bed_line = "\t".join(region_cpt_coords) + "\t".join([gene_name, str(score), region.fields[5]])

        print bed_line

        pybedtools.cleanup(remove_all=True)
        
    if verbose: progress.end()

    #return top_scores(changepoint_scores, 10)
    return

def calc_changepoint(count_vector):
    """Return first changepoint given an IntVector of counts."""
    # Note: Currently only considers single (first) changepoint
    # Can change for multiple (or maybe last?) changepoint

    cpointdata = changepoint.cpt_meanvar(count_vector)
    cpoints_array = changepoint.cpts(cpointdata)

    if len(cpoints_array) != 0:
        cpoint = int(changepoint.cpts(cpointdata)[0])  # parse changepoint values
    else:
        cpoint = 0    # ie, no changepoint found

    return cpoint

def top_scores(score_dict, n):
    """Return top n score entries from score dictionary"""

    score_list = sorted(score_dict.iteritems(), key=operator.itemgetter(1))
    score_list.reverse()
    top_score_list = score_list[:n]

    return top_score_list

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

    pybedtools.set_tempdir('/vol2/home/speach/tmp')

    return changepoints(args.gene_bed, args.signal_bedgraph, args.verbose)

if __name__ == '__main__':
    sys.exit(main())

