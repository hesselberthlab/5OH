import sys
#import ipdb

filename = sys.argv[1]
distance_downstream = 30

for line in open(filename):
    line = line.strip().split("\t")
    chrom, start, stop, gene = line[:4]
    strand = line[10]

    if strand == "+":
        gene_start = int(line[6])
        codon_position = (int(start) - gene_start)%3
        start = int(start) - codon_position
        stop = int(stop) + distance_downstream - codon_position
    
    if strand == "-":
        gene_start = int(line[7])
        codon_position = (int(stop) - gene_start)%3
        stop = int(stop) - codon_position
        start = int(start) - distance_downstream - codon_position + 1

    print "\t".join([chrom, str(start), str(stop), gene+strand])
