#! /usr/bin/env python

'''dna_to_aa: convert DNA sequence to amino acids;
    mostly from: http://www.petercollingridge.co.uk/book/export/html/474'''

import sys

__version__ = '0.1'
__author__ = '@speach'

BASES = ['t', 'c', 'a', 'g']
CODONS = [a+b+c for a in BASES for b in BASES for c in BASES]
AMINO_ACIDS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON_TABLE = dict(zip(CODONS, AMINO_ACIDS))

def translate_tabbed_fasta(filename):
    for line in open(filename):
        gene, dna_sequence = line.strip().split("\t")
        protein_sequence = translate(dna_sequence)
        
        print "\t".join([gene, protein_sequence])

def translate(seq):
    seq = seq.lower().replace('\n', '').replace(' ', '')
    peptide = ''
    
    for i in xrange(0, len(seq), 3):
        codon = seq[i: i+3]
        amino_acid = CODON_TABLE.get(codon, '*')
        if amino_acid != '*':
            peptide += amino_acid
        else:
            break
                
    return peptide


def main():

    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(description=__doc__,
                            version=__version__,
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--filename', help='input filename')

    args = parser.parse_args()

    return translate_tabbed_fasta(args.filename)

if __name__ == '__main__':
    sys.exit(main())
