#! /usr/bin/env python

'''
feature_density: calculate density of signals relative to features.


https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#
bp3-plot-transcription-factor-occupancy-surrounding-the-transcription-start-site

'''
import ipdb
import sys
from itertools import groupby

from toolshed import reader
from pybedtools import BedTool, Interval

__version__ = '$Revision$'

BEDGRAPH_STRANDS = ('+', '-')
SIGNAL_COLNUM = 5
GROUP_COLNUM = 5
SUMMARY_COLNUM = 7

def feature_density(feature_bed_filename,
                    chromsize_filename,
                    pos_bedgraph_filename,
                    neg_bedgraph_filename,
                    invert_strand, flank_size, feature_label,
                    sample_label, library_type,
                    window_resolution, map_operation, group_operation,
                    verbose):

    feature_bedtool = BedTool(feature_bed_filename)
    validate_features(feature_bedtool)

    pos_signal_bedtool = BedTool(pos_bedgraph_filename)
    neg_signal_bedtool = BedTool(neg_bedgraph_filename)

    signal_bedtool = add_strand_to_bedgraph(pos_signal_bedtool,
                                            neg_signal_bedtool,
                                            verbose)

    feature_slop = feature_bedtool.slop(b=flank_size,
                                        g=chromsize_filename)

    feature_windows = make_windows(feature_slop, window_resolution,
                                   invert_strand, verbose)

    feature_map = make_map(feature_windows,
                           signal_bedtool,
                           map_operation,
                           invert_strand,
                           verbose)

    if verbose:
        print >>sys.stderr, ">> grouping data ..."

    feature_grouped = feature_map.groupby(g=GROUP_COLNUM,
                                          c=SUMMARY_COLNUM,
                                          o=group_operation)

    write_table(feature_grouped, flank_size, window_resolution,
                invert_strand, feature_label, sample_label, library_type,
                verbose)

def validate_features(bedtool):
    ''' validate feature BED file.

        Returns:
            Nothing
    ''' 
    named_features = True
    for row in bedtool:
        # name is 4th field
        if row.fields[3] == '.':
            named_features = False

    if not named_features:
        print >>sys.stderr, ">> error: unnamed features in BED file"
        sys.exit(1)

def add_strand_to_bedgraph(pos_bedtool, neg_bedtool, verbose):

    ''' rewrites bedgraph to bed6 format, adding strand.
    
        Returns:
            BedTool
    '''

    if verbose:
        print >>sys.stderr, ">> adding strand to bedgraph data ..."

    groups = ((pos_bedtool, '+'),
              (neg_bedtool, '-'))
    
    intervals = []

    for bedtool, strand in groups:
        for row in bedtool:
            chrom, start, stop, count = row.fields
            fields = [chrom, int(start), int(stop), '.', count, strand]
            intervals.append(Interval(*fields))

    return BedTool(intervals)

def write_table(grouped_bedtool, flank_size, window_resolution,
                invert_strand, feature_label,
                sample_label, library_type, verbose):
    '''Print results in tabular format

        Returns:
            Nothing
    '''
    header_fields = ('#pos', 'rel.pos', 'signal',
                     'library.type', 'feature.label',
                     'sample.label')
    print '\t'.join(header_fields)

    # load the data
    fname = grouped_bedtool.TEMPFILES[-1]
    data = []

    for row in reader(fname, header=['pos','signal']):
        fields = (row['pos'], row['signal'])
        data.append(fields)

    xscale = make_x_scale(flank_size, window_resolution)

    for relpos, (pos, signal) in zip(xscale, data):
        
        fields = [pos, relpos, signal, library_type, feature_label,
                  sample_label]
        print '\t'.join(map(str, fields))

def make_x_scale(flank_size, window_resolution):
    ''' make the x scale for relative positions. for flank_size 500 and
        window_resolution, scale is: -500, -490 ... 0, ... 490, 500

        Returns:
            Iterable
    '''
    xstart = -1 * flank_size
    xend = flank_size + window_resolution
    xscale = range(xstart, xend, window_resolution)

    return xscale

def make_map(windows_bedtool, signal_bedtool, map_operation,
             invert_strand, verbose):
    ''' maps signals onto bed regions using map_operation. mapping takes
        strand comparison into account.

        Returns:
            BedTool, sorted by window number
    '''
    if verbose:
        print >>sys.stderr, ">> making map ... "

    args = {'b':signal_bedtool,
            'c':SIGNAL_COLNUM,
            'o':map_operation,
            'null':0}

    # -s: Require same strandedness
    # -S: Require different strandedness

    if invert_strand:
        args.update({'S':True})
    else:
        args.update({'s':True})

    feature_map = windows_bedtool.map(**args)

    # sort by the window number
    def keyfunc(interval):
        return int(interval.fields[GROUP_COLNUM-1])

    return BedTool(sorted(feature_map, key=keyfunc))

def make_windows(bedtool, window_resolution, invert_strand, verbose):
    ''' 
    Make windows based on specified window_resultion
    Returns:
        BedTool
    '''
    if verbose:
        print >>sys.stderr, ">> making windows ... "

    # select regions from bedtool that have matching strand. have to do
    # this here because make_windows() drops the strand field

    intervals = []
    window_intervals = []
    for strand in BEDGRAPH_STRANDS:
        for interval in bedtool:

            if not invert_strand and interval.strand == strand:
                intervals.append(interval)

            elif invert_strand and interval.strand != strand:
                intervals.append(interval)

        stranded_bedtool = BedTool(intervals)

        # generate the initial set of windows
        args = {'b':stranded_bedtool,
                'w':window_resolution,
                'i':'srcwinnum'}

        windows = BedTool().window_maker(**args)

        # split the 4th field to extract the original name
        # and get the window number, and add the strand for the feature
        for window in windows:
            chrom, start, end, win_field = window.fields

            name, winnum = win_field.split('_')
            fields = [chrom, int(start), int(end),
                      name, winnum, strand]

            window_intervals.append(Interval(*fields))

    def groupkeyfunc(interval):
        # group by name and strand
        return tuple([interval.fields[3], interval.fields[5]])

    grouped_intervals = groupby(window_intervals, key=groupkeyfunc)

    reorder_intervals = []

    for group, intervals in grouped_intervals:

        name, strand = group

        if strand == '+':
            reorder_intervals.extend(intervals)

        elif strand == '-':
            rev_intervals = reverse_interval_windows(intervals)
            reorder_intervals.extend(rev_intervals)

    def sortkeyfunc(interval):
        return tuple([interval.chrom, interval.start])

    return BedTool(sorted(reorder_intervals, key=sortkeyfunc))

def reverse_interval_windows(intervals):
    '''
    Returns:
        list of Intervals
    '''
    # 5th field is the window number
    intervals = tuple(intervals)

    rev_winnums = reversed([interval.fields[4] for interval in intervals])

    reordered = []
    for interval, winnum in zip(intervals, rev_winnums):
        interval_name = interval.fields[3]
        new_interval = Interval(interval.chrom, interval.start,
                                interval.end, interval_name,
                                winnum, interval.strand)

        reordered.append(new_interval)

    return reordered

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... FEATURE_BED SIGNAL_BEDGRAPH CHROM_SIZE"
    version = "%%prog %s" % __version__
    description = ("")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Required")

    group.add_option("-p", "--pos-bedgraph", action="store", type='str',
        default=None, help="pos strand bedgraph signal "
                           "(default: %default)")

    group.add_option("-n", "--neg-bedgraph", action="store", type='str',
        default=None, help="neg strand bedgraph signal "
                           "(default: %default)")

    parser.add_option_group(group)

    group = OptionGroup(parser, "Variables")

    group.add_option("--invert-strand", action="store_true", 
        default=False, help="invert strand signal match if mismatched " 
                           "(default: %default)")

    group.add_option("--flank-size", action="store", type='int',
        default=500, help="flank size (default: %default)")

    group.add_option("--window-resolution", action="store", type='int',
        default=10, help="window resolution (default: %default)")

    group.add_option("--feature-label", action="store", type='str',
        default='feature.label', help="feature label, reported in output"
        " (default: %default)")

    group.add_option("--sample-label", action="store", type='str',
        default='sample.label', help="sample label, reported in output"
        " (default: %default)")

    group.add_option("--library-type", action="store", type='str',
        default='libtype', help="library type, reported in output"
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

    if not options.pos_bedgraph or not options.neg_bedgraph:
        parser.error("specify stranded bedgraph data")

    if len(args) != 2:
        parser.error("specify 32required files")

    if options.flank_size % options.window_resolution != 0:
        parser.error("uneven window spacing, adjust params")

    return options, args

def main(args=sys.argv[1:]):

    options, args = parse_options(args)

    kwargs = {'pos_bedgraph_filename':options.pos_bedgraph,
              'neg_bedgraph_filename':options.neg_bedgraph,
              'invert_strand':options.invert_strand,
              'flank_size':options.flank_size,
              'window_resolution':options.window_resolution,
              'feature_label':options.feature_label,
              'sample_label':options.sample_label,
              'library_type':options.library_type,
              'map_operation':options.map_operation,
              'group_operation':options.group_operation,
              'verbose':options.verbose}

    feature_bed_filename = args[0]
    chromsize_filename = args[1]

    return feature_density(feature_bed_filename, 
                           chromsize_filename,
                           **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
