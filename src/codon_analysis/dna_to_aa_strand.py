#! /usr/bin/env python

'''dna_to_aa: convert DNA sequence to amino acids;
    mostly from: http://www.petercollingridge.co.uk/book/export/html/474'''
    
# Currently cuts off one aa of negative strand data with 3rd codon info

import sys
#import ipdb

__version__ = '0.1'
__author__ = '@speach'

BASES = ['t', 'c', 'a', 'g']
COMP_BASES = ['a', 'g', 't', 'c']
CODONS = [a+b+c for a in BASES for b in BASES for c in BASES]
ANTICODONS = [c+b+a for a in COMP_BASES for b in COMP_BASES for c in COMP_BASES]
AMINO_ACIDS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON_TABLE = dict(zip(CODONS, AMINO_ACIDS))
ANTICODON_TABLE = dict(zip(ANTICODONS, AMINO_ACIDS))

def translate_tabbed_fasta(filename):
    #print "\t".join(["gene","aa1","aa2","aa3","aa4","aa5","aa6","aa7","aa8","aa9","aa10"])
    for line in open(filename):
        gene_and_strand, dna_sequence = line.strip().split("\t")
        gene = gene_and_strand[:-1]
        strand = gene_and_strand[-1]
        protein_sequence = translate(dna_sequence, strand)
        
        print "\t".join([gene, protein_sequence])

def translate(seq, strand):
    # parse seq, initiate variables
    seq = seq.lower().replace('\n', '').replace(' ', '')
    peptide = ''
    iteration = xrange(0, len(seq), 3)
    
    # if positive, use standard codon table
    if strand == "+":
        table = CODON_TABLE
    
    # if negative, use anticodon table and reverse iteration direction
    if strand == "-":
        table = ANTICODON_TABLE
        iteration = reversed(iteration)
    
    for i in iteration:
        codon = seq[i: i+3]
        amino_acid = table.get(codon, 'x')
        if amino_acid == '*':
            peptide += amino_acid
            break
        elif amino_acid != 'x':
            peptide += amino_acid
        else:
            break
        peptide += "\t"
              
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
