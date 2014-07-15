#! /usr/bin/env python

'''
feature_density: calculate density of signals relative to features.


https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#
bp3-plot-transcription-factor-occupancy-surrounding-the-transcription-start-site

'''
import ipdb
import sys

from pybedtools import BedTool, Interval

__version__ = '$Revision$'

BEDGRAPH_SIGNAL_COLNUM = 4
GROUP_COLNUM = 5
SUMMARY_COLNUM = 7

def feature_density(feature_bed_filename, signal_bedgraph_filename,
                    chromsize_filename, window_size, feature_label, 
                    window_resolution, map_operation, group_operation,
                    verbose):

    feature_bedtool = BedTool(feature_bed_filename)
    signal_bedtool = BedTool(signal_bedgraph_filename)

    feature_slop = feature_bedtool.slop(b=window_size,
                                        g=chromsize_filename)

    feature_windows = make_windows(feature_slop, window_resolution, verbose)

    feature_map = make_map(feature_windows, signal_bedtool, map_operation,
                           verbose)

    feature_grouped = feature_map.groupby(g=GROUP_COLNUM,
                                          c=SUMMARY_COLNUM,
                                          o=group_operation)

    pdb.set_trace()

    for row in feature_grouped:
        print '\t'.join(row)

def make_map(windows_bedtool, signal_bedtool, map_operation, verbose):

    if verbose:
        print >>sys.stderr, ">> making map ... "

    feature_map = windows_bedtool.map(b=signal_bedtool,
                                      c=BEDGRAPH_SIGNAL_COLNUM,
                                      o=map_operation,
                                      null=0)
  
    def keyfunc(interval):
        return int(interval.fields[GROUP_COLNUM-1])

    return BedTool(sorted(feature_map, key=keyfunc))

def make_windows(bedtool, window_resolution, verbose):
    ''' 
        Returns:
            BedTool
    '''
    if verbose:
        print >>sys.stderr, ">> making windows ... "

    # generate the initial set of windows
    windows = BedTool().window_maker(b=bedtool, w=window_resolution,
                                     i='srcwinnum')

    # sort by chrom, start and split the 4th field by "_"
    results = []
    for window in windows:
        fs = window.fields
        name, winnum = fs[3].split('_')
        fields = [fs[0], int(fs[1]), int(fs[2]), name, winnum]
        results.append(fields)

    def keyfunc(fs):
        return tuple([fs[0], fs[1]])

    intervals = [Interval(*i) for i in sorted(results, key=keyfunc)]

    return BedTool(intervals)

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... FEATURE_BED SIGNAL_BEDGRAPH CHROM_SIZE"
    version = "%%prog %s" % __version__
    description = ("")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("--window-size", action="store", type='int',
        default=1000, help="window size (default: %default)")

    group.add_option("--window-resolution", action="store", type='int',
        default=10, help="window resolution (default: %default)")

    group.add_option("--feature-label", action="store", type='str',
        default='label', help="feature label, reported in output"
        " (default: %default)")

    group.add_option("--map-operation", action="store", type='str',
        default='mean', help="map operation"
        " (default: %default)")

    group.add_option("--group-operation", action="store", type='str',
        default='sum', help="group operation"
        " (default: %default)")

    group.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if len(args) != 3:
        parser.error("specify 3 required files")

    return options, args

def main(args=sys.argv[1:]):

    options, args = parse_options(args)

    kwargs = {'window_size':options.window_size,
              'window_resolution':options.window_resolution,
              'feature_label':options.feature_label,
              'map_operation':options.map_operation,
              'group_operation':options.group_operation,
              'verbose':options.verbose}

    feature_bed_filename = args[0]
    signal_bedgraph_filename = args[1]
    chromsize_filename = args[2]

    return feature_density(feature_bed_filename, signal_bedgraph_filename,
                           chromsize_filename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
