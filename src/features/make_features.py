#! /usr/bin/env python

'''
make_features: generate feature types from a BED file of genes 


'''
import ipdb
import sys

from pybedtools import BedTool, Interval

BED_FIELDNAMES = ('chrom', 'start', 'end', 'name', 'score',
                  'strand', 'thickStart', 'thickEnd',
                  'itemRgb', 'blockCount', 'blockSizes',
                  'blockStarts')

__version__ = '$Revision$'

def make_features(feature_bed_filename, verbose):

    feature_bedtool = BedTool(feature_bed_filename)

def make_transcription_starts(bedtool):
    ''' '''
    intervals = []
    for interval in bedtool:
        if interval.strand == '+':
        else:
    return BedTool(intervals)

def make_transcription_stops(bedtool):
    intervals = []
    for interval in bedtool:
        if interval.strand == '+':
        else:
    return BedTool(intervals)

def make_coding_starts(bedtool):
    intervals = []
    for interval in bedtool:
        if interval.strand == '+':
        else:
    return BedTool(intervals)

def make_coding_stops(bedtool):
    intervals = []
    for interval in bedtool:
        if interval.strand == '+':
        else:
    return BedTool(intervals)

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... FEATURE_BED "
    version = "%%prog %s" % __version__
    description = ("")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if len(args) != 1:
        parser.error("specify BED file")

    return options, args

def main(args=sys.argv[1:]):

    options, args = parse_options(args)

    kwargs = {'verbose':options.verbose}

    feature_bed_filename = args[0]

    return make_features(feature_bed_filename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
