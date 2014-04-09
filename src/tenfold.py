#! /usr/bin/env python

import RNA
import re
import sys

def tenfold(fasta):
    """Fold tabbed FASTA input (gene;chr:start-stop \t sequence) and
    convert to binary output representing highest likelihood secondary  
    structure
    1 = double-stranded
    0 = single-standed
    """
    for line in open(fasta):
        label, sequence = line.strip().split("\t")
        dotplot, fe = RNA.fold(sequence)
        re1 = re.sub(r'\(|\)',r'1',dotplot)
        re10 = re.sub(r'\.',r'0',re1)
        print label + "\t" + re10

tenfold(sys.argv[1])
