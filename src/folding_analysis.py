#! /usr/bin/env python

'''5OH-folding-analysis: XXX
'''
import sys
from itertools import izip

from pybedtools import BedTool
from RNA import fold, pf_fold

def folding_analysis(bedfilename, fastafilename, verbose):

    bedtool = BedTool(bedfilename)
    for region in bedtool:

        region_seq = bedtool.sequence()

        struct, mfe = RNA.fold(region_seq) 

        for pos, nuc  in enumerate(region_seq):
            struct_char = struct[pos]

def is_paired(char):
    '''
    '''
    if char in '()':
        return True
    return False

def is_sstranded(char):
    '''
    '''
    if char in '.':
        return True
    return False

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... BEDFILENAME FASTAFILENAME"
    version = "%%prog %s" % __version__
    description = ("")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    options, args = parser.parse_args(args)

    if len(args) < 1:
        parser.error("specify genomedata filename")
   
    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)
    bedfilename = args[0]
    fastafilename = args[1]
    kwargs = {'verbose':options.verbose}
    return folding_analysis(bedfilename, fastafilename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
