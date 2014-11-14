'''Input: fasta file generated from sacCer1 w/ poly-aa_genes.bed 100bp stretches
   Output: tab-delimited gene, motif, and codons'''

import sys
from Bio.Seq import Seq # package with useful methods for sequences

filename = sys.argv[1]

gene_seq = []
linecount = 1

def print_gene_seq(gene_seq):
    gene, motif, strand, seq = gene_seq
    
    if strand == "+":
        seq = Seq(seq[50:])
        
    if strand == "-":
        seq = Seq(seq[:51]).reverse_complement()
    
    # Parsing motif information
    if motif[-2:] == '15':
        motif_len = motif[-2:]
    elif motif[-1] == '5':
        motif_len = 5
    else:
        motif_len = motif[-2:]
    
    codons = []
    #for i in range(0,int(motif_len)): # full call to include all calculated poly-aas
    for i in range(0,5): # just including the poly-aas of length 5.
        codon = str(seq[i*3:i*3+3]) # need to cast Seq object as string
        codons.append(codon)
    
    print "\t".join([gene, motif, "\t".join(codons)])
    
def invert_seq(seq):
    seq = seq[::-1] # reverse string

print "\t".join(["gene", "motif","codon.1","codon.2","codon.3","codon.4","codon.5"])

for line in open(filename):
    line = line.strip()

    if linecount % 2 == 1: # if name line
        name, motif, strand = line.lstrip(">").split(":")
        gene_seq.append(name)
        gene_seq.append(motif)
        gene_seq.append(strand)
    
    else: # line is seq
        seq = line
        gene_seq.append(seq)
        print_gene_seq(gene_seq)

        gene_seq = []

    linecount += 1



    
