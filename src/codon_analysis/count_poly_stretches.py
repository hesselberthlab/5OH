#! /usr/bin/env python

'''count_poly_stretches.py: count number of genes in sarCer with poly-amino acid stretches'''

from subprocess import Popen, PIPE

aminos = ["D", "E", "[DE]", "R", "K", "[RK]"]
protein_list = '/vol2/home/speach/ref/genomes/sacCer1/orf_trans.tab'

def count_poly_stretches(aminos, protein_list, mismatch=0):
    # set shell command call; agrep allows for a regex with a predifined number of mismatches
    if mismatch == 0:
        grep = "grep"
    else:
        grep = "agrep -%s" % mismatch
    
    # loop through amino acids, grep genes via bash subprocess
    for i in range (1,5):
        for amino in aminos:
            polyaa = amino * 5 * i
            
            bashcommand = "%s %s %s | wc" % (grep, polyaa, protein_list)
            process = Popen(bashcommand.split(), stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()
            num_genes = len(stdout.split("\n")) - 1
        
            print "\t".join([str(mismatch), str(num_genes), amino, str(5*i)])

print "\t".join(["mismatch", "num.genes", "polyaa", "num.rep"])
for i in range(0,5):
    count_poly_stretches(aminos, protein_list, i)
    
