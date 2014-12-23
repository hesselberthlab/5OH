#! /usr/bin/env python

for line in open('codon_usage_yeast.tab'):
    line = line.strip().split('\t')
    codon_data_len = len(line)/4
    for i in range(1,5):
        codon_data = line[(i-1)*4:i*4+1]
        print codon_data
        print "\t".join([codon_data[0], codon_data[2]])
    
