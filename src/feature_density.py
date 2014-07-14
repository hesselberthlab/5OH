#! /usr/bin/env python

'''
feature_density: calculate density of signals relative to features.


https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#
bp3-plot-transcription-factor-occupancy-surrounding-the-transcription-start-site

'''
import sys
from pybedtools import BedTool

__version__ = '$Revision$'

def feature_density(feature_bedfilename, signal_bedgraph_filename,
                    chromsize_filename, flank_size, window_resolution,
                    map_operation, group_operation, verbose):

    feature_bedtool = BedTool(feature_bedfile)
    signal_bedtool = BedTool(signal_bedgraphfile)

    feature_slop = feature_bedtool.slop(b=window_size,
                                        g=chromsize_filename)

    feature_windows = feature_slop.makewindows(w=window_resolution,
                                               i=srcwinnum)

    feature_map = feature_windows.map(c=signal_colnum, o=map_oper, null=0)

    group_col = 5
    summary_col = 6
    feature_grouped = feature_map.groupby(g=group_col, c=summary_col,
                                          o=group_oper)

def make_feature_windows(feature_bedtool, window_resolution, verbose):

    feature_windows = feature_bedtool.makewindows(w=window_resolution,
                                                  i=srcwinnum)
    return result

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... FEATURE_BED SIGNAL_BEDGRAPH, CHROM_SIZE"
    version = "%%prog %s" % __version__
    description = ("")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("--flank-size", action="store", type='int',
        default=500, help="flank size (default: %default)")

    group.add_option("--window-resolution", action="store", type='int',
        default=10, help="window resolution (default: %default)")

    group.add_option("--feature-label", action="store", type='str',
        default='', help="feature label, reported in output"
        " (default: %default)")

    group.add_option("--map-operation", action="store", type='str',
        default='sum', help="map operation"
        " (default: %default)")

    group.add_option("--group-operation", action="store", type='str',
        default='mean', help="group operation"
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

    kwargs = {'flank_size':options.flank_size,
              'window_res':options.window_resolution,
              'feature_label':options.feature_label,
              'map_oper':options.map_oper,
              'group_oper':options.group_oper,
              'verbose':options.verbose}

    feature_bedfile = args[0]
    signal_bedgraphfile = args[1]
    chromsize_fielname = args[3]

    return feature_density(feature_bedfilename, signal_bedgraphfilename,
                           chromsize_filename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
