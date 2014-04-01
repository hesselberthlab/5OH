#! /usr/bin/env python

import sys
from pysam import Samfile

def junction_anlaysis(stub):
    pass

def parse_samfile(samfilename):
    ''' Function to parse UMI and read sequence from sam (bam) file
       
        Params:
            samfilename (str): sam filename
        Returns:
            UMI (str), sequence (str)
    '''
    with Samfile(samfilename, 'rb') as samfile:
        pass 
