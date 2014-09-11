#! /usr/bin/env python

''' changpoints: calc changepoints '''

import sys
import pdb
import rpy2.robjects as robjects

from pybedtools import BedTool
from rpy2.robjects.packages import importr

__version__ = '0.1'

# Import Changepoint package from R
try:
    changepoint = importr('changepoint')
except ImportError:
    print >>sys.stderr, "install R changepoint package"
    sys.exit(1)

def changepoints(gene_bed, signal_bedgraph):

    genes = BedTool(gene_bed)
    signal = BedTool(signal_bedgraph)

    for region in genes:      # loop through each gene
        region_bedtool = BedTool([region.fields[:3]])
        signal_data = signal.intersect(region_bedtool)

        # Create list from count intersection
        count_list = [int(datum.fields[3]) for datum in signal_data]

        # Convert to Int Vector
        count_vector = robjects.IntVector(count_list)

        region_cpt = calc_changepoint(count_vector)

        # Todo: Store as tuple with region entry.

def calc_changepoint(count_vector):
    """Return first changepoint given an IntVector of counts."""
    # Note: Currently only considers single (first) changepoint
    # Can change for multiple (or maybe last?) changepoint

    cpointdata = changepoint.cpt_meanvar(count_vector)
    cpoint = int(changepoint.cpts(cpoint)[0])  # parse changepoint values

    return cpoint

def main():

    from argparse import ArgumentParser, RawDescriptionHelpFormatter

    parser = ArgumentParser(description=__doc__,
                            version=__version__,
                            formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument('gene_bed', help='BED file of genes')
    parser.add_argument('signal_bedgraph', help='bedgraph signal')

    args = parser.parse_args()

    return changepoints(args.gene_bed, args.signal_bedgraph)

if __name__ == '__main__':
    sys.exit(main())

